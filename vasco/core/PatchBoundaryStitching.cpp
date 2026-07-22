#include "vasco/core/PatchBoundaryStitching.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <CGAL/AABB_segment_primitive_3.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>

namespace vasco::patch_boundary_stitching {
namespace {

using SlicerVec3 = vasco::Slicer::Vec3;
using Vector3 = Eigen::Vector3d;
using Tri3 = vasco::Slicer::Tri3;
using EdgeKey = std::pair<int, int>;
using AabbKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using AabbPoint = AabbKernel::Point_3;
using AabbSegment = AabbKernel::Segment_3;
using AabbBox = AabbKernel::Iso_cuboid_3;
using AabbSegmentContainer = std::vector<AabbSegment>;
using AabbSegmentIterator = AabbSegmentContainer::const_iterator;
using AabbPrimitive = CGAL::AABB_segment_primitive_3<AabbKernel, AabbSegmentIterator>;
using AabbTraits = CGAL::AABB_traits_3<AabbKernel, AabbPrimitive>;
using AabbTree = CGAL::AABB_tree<AabbTraits>;
using ValidationMesh = CGAL::Surface_mesh<AabbPoint>;

struct EdgeKeyHash {
	std::size_t operator()(const EdgeKey& edge) const
	{
		const std::size_t first = std::hash<int>{}(edge.first);
		const std::size_t second = std::hash<int>{}(edge.second);
		return first ^ (second + 0x9e3779b9u + (first << 6) + (first >> 2));
	}
};

struct EdgeRecord {
	int first_face_id = -1;
	std::size_t use_count = 0;

	void AddFace(int face_id)
	{
		if (use_count == 0) {
			first_face_id = face_id;
		}
		++use_count;
	}
};

struct BoundaryEdge {
	EdgeKey key;
	int face_id = -1;
	int source_patch_id = -1;
};

struct Topology {
	std::vector<BoundaryEdge> boundary_edges;
	std::vector<std::vector<int>> incident_faces;
	std::vector<std::vector<int>> incident_boundary_edges;
	bool has_non_manifold_edge = false;
};

struct Candidate {
	int point_id = -1;
	int edge_id = -1;
	double distance = 0.0;
	double parameter = 0.0;
	double score = 0.0;
	Vector3 projected_point = Vector3::Zero();
};

Vector3 ToVector3(const SlicerVec3& point)
{
	return { point[0], point[1], point[2] };
}

SlicerVec3 ToSlicerVec3(const Vector3& point)
{
	return { point.x(), point.y(), point.z() };
}

EdgeKey MakeEdgeKey(int a, int b)
{
	return (a < b) ? EdgeKey{ a, b } : EdgeKey{ b, a };
}

// Computes a finite unit normal and rejects invalid indices and degenerate faces.
bool ComputeFaceNormal(const vasco::Slicer& mesh, const Tri3& triangle, Vector3& normal)
{
	if (triangle[0] < 0 || triangle[1] < 0 || triangle[2] < 0
		|| triangle[0] >= static_cast<int>(mesh.positions.size())
		|| triangle[1] >= static_cast<int>(mesh.positions.size())
		|| triangle[2] >= static_cast<int>(mesh.positions.size())) {
		return false;
	}
	const Vector3 point0 = ToVector3(mesh.positions[triangle[0]]);
	const Vector3 point1 = ToVector3(mesh.positions[triangle[1]]);
	const Vector3 point2 = ToVector3(mesh.positions[triangle[2]]);
	if (!point0.allFinite() || !point1.allFinite() || !point2.allFinite()) {
		return false;
	}
	const Vector3 ab = point1 - point0;
	const Vector3 ac = point2 - point0;
	normal = ab.cross(ac);
	const double length = normal.norm();
	if (!std::isfinite(length) || length <= 1e-24) {
		return false;
	}
	normal /= length;
	return true;
}

double ComputeModelDiagonal(const vasco::Slicer& mesh)
{
	if (mesh.positions.empty()) {
		return 0.0;
	}
	Vector3 minimum = ToVector3(mesh.positions.front());
	Vector3 maximum = minimum;
	for (const auto& point : mesh.positions) {
		const Vector3 vector = ToVector3(point);
		if (!vector.allFinite()) {
			return std::numeric_limits<double>::quiet_NaN();
		}
		minimum = minimum.cwiseMin(vector);
		maximum = maximum.cwiseMax(vector);
	}
	return (maximum - minimum).norm();
}

double ComputeModelTolerance(double model_diagonal, const StitchOptions& options)
{
	return std::clamp(
		model_diagonal * options.model_relative_tolerance,
		options.absolute_min_tolerance,
		options.absolute_max_tolerance);
}

double ComputeEdgeTolerance(
	double edge_length,
	double model_tolerance,
	const StitchOptions& options)
{
	return std::clamp(
		edge_length * options.edge_relative_tolerance,
		options.absolute_min_tolerance,
		model_tolerance);
}

bool AreOptionsValid(const StitchOptions& options)
{
	return std::isfinite(options.model_relative_tolerance)
		&& std::isfinite(options.edge_relative_tolerance)
		&& std::isfinite(options.absolute_min_tolerance)
		&& std::isfinite(options.absolute_max_tolerance)
		&& std::isfinite(options.normal_cosine_threshold)
		&& std::isfinite(options.boundary_direction_cosine_threshold)
		&& options.model_relative_tolerance > 0.0
		&& options.edge_relative_tolerance > 0.0
		&& options.absolute_min_tolerance > 0.0
		&& options.absolute_min_tolerance <= options.absolute_max_tolerance
		&& options.normal_cosine_threshold >= 0.0
		&& options.normal_cosine_threshold <= 1.0
		&& options.boundary_direction_cosine_threshold >= 0.0
		&& options.boundary_direction_cosine_threshold <= 1.0
		&& options.max_iterations > 0;
}

// Builds edge incidence and the boundary-only adjacency needed for candidate tests.
bool BuildTopology(
	const vasco::Slicer& mesh,
	const std::vector<int>& face_source_patch_ids,
	Topology& topology)
{
	topology = {};
	if (mesh.triangles.size() != face_source_patch_ids.size()) {
		return false;
	}
	topology.incident_faces.resize(mesh.positions.size());
	std::unordered_map<EdgeKey, EdgeRecord, EdgeKeyHash> edge_uses;
	edge_uses.reserve(mesh.triangles.size() * 2 + 1);
	for (int face_id = 0; face_id < static_cast<int>(mesh.triangles.size()); ++face_id) {
		const auto& triangle = mesh.triangles[face_id];
		Vector3 normal;
		if (!ComputeFaceNormal(mesh, triangle, normal)) {
			return false;
		}
		for (int corner = 0; corner < 3; ++corner) {
			topology.incident_faces[triangle[corner]].push_back(face_id);
			const int from = triangle[corner];
			const int to = triangle[(corner + 1) % 3];
			edge_uses[MakeEdgeKey(from, to)].AddFace(face_id);
		}
	}

	for (const auto& entry : edge_uses) {
		if (entry.second.use_count > 2) {
			topology.has_non_manifold_edge = true;
			continue;
		}
		if (entry.second.use_count != 1) {
			continue;
		}
		topology.boundary_edges.push_back({
			entry.first,
			entry.second.first_face_id,
			face_source_patch_ids[entry.second.first_face_id]
		});
	}
	std::sort(
		topology.boundary_edges.begin(),
		topology.boundary_edges.end(),
		[](const BoundaryEdge& lhs, const BoundaryEdge& rhs) {
			return lhs.key < rhs.key;
		});

	topology.incident_boundary_edges.resize(mesh.positions.size());
	for (int edge_id = 0; edge_id < static_cast<int>(topology.boundary_edges.size()); ++edge_id) {
		const auto& edge = topology.boundary_edges[edge_id];
		topology.incident_boundary_edges[edge.key.first].push_back(edge_id);
		topology.incident_boundary_edges[edge.key.second].push_back(edge_id);
	}
	return true;
}

class BoundaryEdgeAabbTree {
public:
	BoundaryEdgeAabbTree(
		const std::vector<BoundaryEdge>& edges,
		const std::vector<SlicerVec3>& positions)
	{
		segments_.reserve(edges.size());
		for (const auto& edge : edges) {
			const auto& begin = positions[edge.key.first];
			const auto& end = positions[edge.key.second];
			segments_.emplace_back(
				AabbPoint(begin[0], begin[1], begin[2]),
				AabbPoint(end[0], end[1], end[2]));
		}
		if (!segments_.empty()) {
			tree_.insert(segments_.cbegin(), segments_.cend());
			tree_.build();
		}
	}

	void Query(const Vector3& point, double radius, std::vector<int>& edge_ids) const
	{
		edge_ids.clear();
		if (segments_.empty()) {
			return;
		}
		const AabbPoint minimum(
			point.x() - radius,
			point.y() - radius,
			point.z() - radius);
		const AabbPoint maximum(
			point.x() + radius,
			point.y() + radius,
			point.z() + radius);
		primitive_ids_.clear();
		tree_.all_intersected_primitives(
			AabbBox(minimum, maximum),
			std::back_inserter(primitive_ids_));
		edge_ids.reserve(primitive_ids_.size());
		for (const auto primitive_id : primitive_ids_) {
			edge_ids.push_back(static_cast<int>(std::distance(segments_.cbegin(), primitive_id)));
		}
	}

private:
	AabbSegmentContainer segments_;
	AabbTree tree_;
	mutable std::vector<AabbSegmentIterator> primitive_ids_;
};

// Projects onto the supporting segment line and returns both the parameter and distance.
bool ProjectPointOntoSegment(
	const Vector3& point,
	const Vector3& begin,
	const Vector3& end,
	double& parameter,
	Vector3& projected_point,
	double& distance)
{
	const Vector3 edge = end - begin;
	const double edge_length_squared = edge.squaredNorm();
	if (edge_length_squared <= 1e-24) {
		return false;
	}
	parameter = (point - begin).dot(edge) / edge_length_squared;
	projected_point = begin + parameter * edge;
	distance = (point - projected_point).norm();
	return true;
}

bool HasCompatibleBoundaryDirection(
	int point_id,
	const BoundaryEdge& candidate_edge,
	const Topology& topology,
	const vasco::Slicer& mesh,
	double cosine_threshold)
{
	Vector3 candidate_direction = ToVector3(mesh.positions[candidate_edge.key.second])
		- ToVector3(mesh.positions[candidate_edge.key.first]);
	if (candidate_direction.norm() <= 1e-24) {
		return false;
	}
	candidate_direction.normalize();
	for (int edge_id : topology.incident_boundary_edges[point_id]) {
		const BoundaryEdge& incident_edge = topology.boundary_edges[edge_id];
		const int other = (incident_edge.key.first == point_id)
			? incident_edge.key.second
			: incident_edge.key.first;
		Vector3 incident_direction = ToVector3(mesh.positions[other]) - ToVector3(mesh.positions[point_id]);
		if (incident_direction.norm() > 1e-24
			&& std::abs(candidate_direction.dot(incident_direction.normalized())) >= cosine_threshold) {
			return true;
		}
	}
	return false;
}

bool HasCompatibleFaceNormal(
	int point_id,
	const BoundaryEdge& candidate_edge,
	const Topology& topology,
	const vasco::Slicer& mesh,
	double cosine_threshold)
{
	Vector3 edge_face_normal;
	if (!ComputeFaceNormal(mesh, mesh.triangles[candidate_edge.face_id], edge_face_normal)) {
		return false;
	}
	Vector3 average_normal = Vector3::Zero();
	for (int face_id : topology.incident_faces[point_id]) {
		Vector3 normal;
		if (ComputeFaceNormal(mesh, mesh.triangles[face_id], normal)) {
			average_normal += normal;
		}
	}
	if (average_normal.norm() <= 1e-24) {
		return false;
	}
	return edge_face_normal.dot(average_normal.normalized()) >= cosine_threshold;
}

bool ComesFromNewerPatch(
	int point_id,
	const BoundaryEdge& candidate_edge,
	const Topology& topology,
	const std::vector<int>& face_source_patch_ids)
{
	int newest_point_patch_id = std::numeric_limits<int>::min();
	for (int face_id : topology.incident_faces[point_id]) {
		newest_point_patch_id = std::max(newest_point_patch_id, face_source_patch_ids[face_id]);
	}
	return candidate_edge.source_patch_id < newest_point_patch_id;
}

// Uses the CGAL boundary-edge AABB tree to find the best source-aware T-junction
// candidate for each boundary vertex.
std::vector<Candidate> FindBestCandidates(
	const vasco::Slicer& mesh,
	const std::vector<int>& face_source_patch_ids,
	const Topology& topology,
	double model_tolerance,
	const StitchOptions& options)
{
	BoundaryEdgeAabbTree aabb_tree(topology.boundary_edges, mesh.positions);
	std::vector<Candidate> candidates;
	std::vector<int> nearby_edges;
	for (int point_id = 0; point_id < static_cast<int>(mesh.positions.size()); ++point_id) {
		if (topology.incident_boundary_edges[point_id].empty()) {
			continue;
		}
		aabb_tree.Query(ToVector3(mesh.positions[point_id]), model_tolerance, nearby_edges);
		Candidate best;
		best.score = std::numeric_limits<double>::infinity();
		for (int edge_id : nearby_edges) {
			const BoundaryEdge& edge = topology.boundary_edges[edge_id];
			if (point_id == edge.key.first || point_id == edge.key.second
				|| !ComesFromNewerPatch(point_id, edge, topology, face_source_patch_ids)) {
				continue;
			}
			const double edge_length = (
				ToVector3(mesh.positions[edge.key.second])
				- ToVector3(mesh.positions[edge.key.first])).norm();
			const double tolerance = ComputeEdgeTolerance(edge_length, model_tolerance, options);
			double parameter = 0.0;
			double distance = 0.0;
			Vector3 projected_point;
			if (!ProjectPointOntoSegment(
				ToVector3(mesh.positions[point_id]),
				ToVector3(mesh.positions[edge.key.first]),
				ToVector3(mesh.positions[edge.key.second]),
				parameter,
				projected_point,
				distance)
				|| distance > tolerance
				|| parameter * edge_length <= tolerance
				|| (1.0 - parameter) * edge_length <= tolerance
				|| !HasCompatibleBoundaryDirection(
					point_id,
					edge,
					topology,
					mesh,
					options.boundary_direction_cosine_threshold)
				|| !HasCompatibleFaceNormal(
					point_id,
					edge,
					topology,
					mesh,
					options.normal_cosine_threshold)) {
				continue;
			}

			const double score = distance / tolerance;
			if (score < best.score) {
				best = { point_id, edge_id, distance, parameter, score, projected_point };
			}
		}
		if (best.edge_id >= 0) {
			candidates.push_back(best);
		}
	}
	return candidates;
}

bool TriangleKeepsOrientation(
	const vasco::Slicer& original_mesh,
	const vasco::Slicer& proposed_mesh,
	const Tri3& original_triangle,
	const Tri3& proposed_triangle,
	double normal_cosine_threshold)
{
	Vector3 original_normal;
	Vector3 proposed_normal;
	return ComputeFaceNormal(original_mesh, original_triangle, original_normal)
		&& ComputeFaceNormal(proposed_mesh, proposed_triangle, proposed_normal)
		&& original_normal.dot(proposed_normal) >= normal_cosine_threshold;
}

bool ValidateMovedPointFaces(
	const vasco::Slicer& original_mesh,
	const vasco::Slicer& proposed_mesh,
	const Topology& topology,
	const std::set<int>& moved_points,
	double normal_cosine_threshold)
{
	std::set<int> affected_faces;
	for (int point_id : moved_points) {
		affected_faces.insert(
			topology.incident_faces[point_id].begin(),
			topology.incident_faces[point_id].end());
	}
	for (int face_id : affected_faces) {
		if (!TriangleKeepsOrientation(
			original_mesh,
			proposed_mesh,
			original_mesh.triangles[face_id],
			proposed_mesh.triangles[face_id],
			normal_cosine_threshold)) {
			return false;
		}
	}
	return true;
}

// Splits one triangle along its original oriented boundary edge. The generated
// triangle fan follows the original cyclic vertex order and therefore its winding.
bool BuildSplitTriangles(
	const vasco::Slicer& original_mesh,
	const vasco::Slicer& proposed_mesh,
	const BoundaryEdge& edge,
	std::vector<Candidate> candidates,
	std::vector<Tri3>& split_triangles,
	double normal_cosine_threshold)
{
	const Tri3& original_triangle = original_mesh.triangles[edge.face_id];
	int edge_corner = -1;
	for (int corner = 0; corner < 3; ++corner) {
		if (MakeEdgeKey(
			original_triangle[corner],
			original_triangle[(corner + 1) % 3]) == edge.key) {
			edge_corner = corner;
			break;
		}
	}
	if (edge_corner < 0) {
		return false;
	}

	const int begin = original_triangle[edge_corner];
	const int end = original_triangle[(edge_corner + 1) % 3];
	const int opposite = original_triangle[(edge_corner + 2) % 3];
	const Vector3 direction = ToVector3(proposed_mesh.positions[end])
		- ToVector3(proposed_mesh.positions[begin]);
	const double length_squared = direction.squaredNorm();
	for (auto& candidate : candidates) {
		candidate.parameter = (
			ToVector3(proposed_mesh.positions[candidate.point_id])
			- ToVector3(proposed_mesh.positions[begin])).dot(direction) / length_squared;
	}
	std::sort(candidates.begin(), candidates.end(), [](const Candidate& lhs, const Candidate& rhs) {
		return lhs.parameter < rhs.parameter;
	});

	std::vector<int> edge_vertices;
	edge_vertices.reserve(candidates.size() + 2);
	edge_vertices.push_back(begin);
	for (const auto& candidate : candidates) {
		if (edge_vertices.back() != candidate.point_id) {
			edge_vertices.push_back(candidate.point_id);
		}
	}
	edge_vertices.push_back(end);

	split_triangles.clear();
	for (std::size_t index = 0; index + 1 < edge_vertices.size(); ++index) {
		Tri3 triangle = { edge_vertices[index], edge_vertices[index + 1], opposite };
		if (!TriangleKeepsOrientation(
			original_mesh,
			proposed_mesh,
			original_triangle,
			triangle,
			normal_cosine_threshold)) {
			return false;
		}
		split_triangles.push_back(triangle);
	}
	return split_triangles.size() >= 2;
}

// Converts to CGAL::Surface_mesh for the final structural and manifoldness checks.
bool ValidateFinalMesh(const vasco::Slicer& mesh, const std::vector<int>& face_source_patch_ids)
{
	if (mesh.triangles.size() != face_source_patch_ids.size()) {
		return false;
	}
	ValidationMesh validation_mesh;
	std::vector<ValidationMesh::Vertex_index> vertices;
	vertices.reserve(mesh.positions.size());
	for (const auto& point : mesh.positions) {
		vertices.push_back(validation_mesh.add_vertex(AabbPoint(point[0], point[1], point[2])));
	}

	for (const auto& triangle : mesh.triangles) {
		Vector3 normal;
		if (!ComputeFaceNormal(mesh, triangle, normal)) {
			return false;
		}
		if (validation_mesh.add_face(
			vertices[triangle[0]],
			vertices[triangle[1]],
			vertices[triangle[2]]) == ValidationMesh::null_face()) {
			return false;
		}
	}
	if (!CGAL::is_valid_polygon_mesh(validation_mesh)) {
		return false;
	}
	for (const auto vertex : validation_mesh.vertices()) {
		if (CGAL::Polygon_mesh_processing::is_non_manifold_vertex(vertex, validation_mesh)) {
			return false;
		}
	}
	return validation_mesh.number_of_faces() == mesh.triangles.size();
}

} // namespace

bool StitchPatchBoundariesWithTolerance(
	vasco::Slicer& mesh,
	std::vector<int>& face_source_patch_ids,
	const StitchOptions& options,
	StitchStats& stats)
{
	stats = {};
	if (mesh.triangles.empty()) {
		return true;
	}
	if (!AreOptionsValid(options)
		|| mesh.triangles.size() != face_source_patch_ids.size()) {
		return false;
	}

	const vasco::Slicer original_mesh = mesh;
	const std::vector<int> original_face_sources = face_source_patch_ids;
	stats.model_diagonal = ComputeModelDiagonal(mesh);
	stats.model_tolerance = ComputeModelTolerance(stats.model_diagonal, options);
	if (!std::isfinite(stats.model_diagonal)
		|| !std::isfinite(stats.model_tolerance)
		|| stats.model_diagonal <= 0.0
		|| stats.model_tolerance <= 0.0) {
		return false;
	}

	for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
		Topology topology;
		if (!BuildTopology(mesh, face_source_patch_ids, topology) || topology.has_non_manifold_edge) {
			mesh = original_mesh;
			face_source_patch_ids = original_face_sources;
			return false;
		}
		if (iteration == 0) {
			stats.initial_boundary_edge_count = topology.boundary_edges.size();
		}

		const std::vector<Candidate> candidates = FindBestCandidates(
			mesh,
			face_source_patch_ids,
			topology,
			stats.model_tolerance,
			options);
		if (iteration == 0) {
			stats.initial_t_junction_count = candidates.size();
		}
		if (candidates.empty()) {
			break;
		}

		std::map<int, std::vector<Candidate>> candidates_by_edge;
		for (const auto& candidate : candidates) {
			candidates_by_edge[candidate.edge_id].push_back(candidate);
		}
		std::map<int, int> selected_edge_by_face;
		for (const auto& entry : candidates_by_edge) {
			const int face_id = topology.boundary_edges[entry.first].face_id;
			auto selected = selected_edge_by_face.find(face_id);
			if (selected == selected_edge_by_face.end()
				|| entry.second.size() > candidates_by_edge[selected->second].size()) {
				selected_edge_by_face[face_id] = entry.first;
			}
		}

		vasco::Slicer proposed_mesh = mesh;
		std::set<int> moved_points;
		std::map<int, std::vector<Tri3>> replacements;
		for (const auto& selected : selected_edge_by_face) {
			const int edge_id = selected.second;
			const auto& edge_candidates = candidates_by_edge[edge_id];
			for (const auto& candidate : edge_candidates) {
				proposed_mesh.positions[candidate.point_id] = ToSlicerVec3(candidate.projected_point);
				moved_points.insert(candidate.point_id);
			}

			std::vector<Tri3> split_triangles;
			if (!BuildSplitTriangles(
				mesh,
				proposed_mesh,
				topology.boundary_edges[edge_id],
				edge_candidates,
				split_triangles,
				options.normal_cosine_threshold)) {
				for (const auto& candidate : edge_candidates) {
					proposed_mesh.positions[candidate.point_id] = mesh.positions[candidate.point_id];
					moved_points.erase(candidate.point_id);
				}
				stats.rejected_candidate_count += edge_candidates.size();
				continue;
			}
			replacements[selected.first] = std::move(split_triangles);
		}

		if (replacements.empty()
			|| !ValidateMovedPointFaces(
				mesh,
				proposed_mesh,
				topology,
				moved_points,
				options.normal_cosine_threshold)) {
			if (!replacements.empty()) {
				stats.rejected_candidate_count += moved_points.size();
			}
			break;
		}

		std::vector<Tri3> new_triangles;
		std::vector<int> new_face_sources;
		new_triangles.reserve(mesh.triangles.size() + candidates.size());
		new_face_sources.reserve(face_source_patch_ids.size() + candidates.size());
		for (int face_id = 0; face_id < static_cast<int>(mesh.triangles.size()); ++face_id) {
			const auto replacement = replacements.find(face_id);
			if (replacement == replacements.end()) {
				new_triangles.push_back(mesh.triangles[face_id]);
				new_face_sources.push_back(face_source_patch_ids[face_id]);
				continue;
			}
			for (const auto& triangle : replacement->second) {
				new_triangles.push_back(triangle);
				new_face_sources.push_back(face_source_patch_ids[face_id]);
			}
			++stats.split_face_count;
			stats.added_face_count += replacement->second.size() - 1;
		}

		mesh.positions.swap(proposed_mesh.positions);
		mesh.triangles.swap(new_triangles);
		face_source_patch_ids.swap(new_face_sources);
		stats.projected_vertex_count += moved_points.size();
		stats.split_edge_count += replacements.size();
		++stats.iteration_count;
	}

	Topology final_topology;
	if (!BuildTopology(mesh, face_source_patch_ids, final_topology)
		|| final_topology.has_non_manifold_edge
		|| !ValidateFinalMesh(mesh, face_source_patch_ids)) {
		mesh = original_mesh;
		face_source_patch_ids = original_face_sources;
		return false;
	}
	stats.final_boundary_edge_count = final_topology.boundary_edges.size();
	stats.remaining_t_junction_count = FindBestCandidates(
		mesh,
		face_source_patch_ids,
		final_topology,
		stats.model_tolerance,
		options).size();
	return true;
}

} // namespace vasco::patch_boundary_stitching
