#include "vasco/core/MeshValidation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>

#include <CGAL/IO/STL.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/boost/graph/helpers.h>

#include "katana/katana.h"

namespace vasco::mesh_validation {
namespace {

namespace PMP = CGAL::Polygon_mesh_processing;

struct PolygonSoupStats {
	std::size_t point_count = 0;
	std::size_t triangle_count = 0;
	std::size_t non_triangle_polygon_count = 0;
	std::size_t invalid_index_polygon_count = 0;
	std::size_t invalid_point_count = 0;
	std::size_t degenerate_triangle_count = 0;
	std::size_t duplicate_triangle_count = 0;
	bool is_polygon_mesh = false;
};

using QuantizedCoordinateKey = std::tuple<long long, long long, long long>;
using TriangleCoordinateKey = std::array<QuantizedCoordinateKey, 3>;

QuantizedCoordinateKey MakeQuantizedCoordinateKeyForSoup(const Point_3& p, double scale = 1e9)
{
	return {
		static_cast<long long>(std::llround(CGAL::to_double(p.x()) * scale)),
		static_cast<long long>(std::llround(CGAL::to_double(p.y()) * scale)),
		static_cast<long long>(std::llround(CGAL::to_double(p.z()) * scale))
	};
}

bool HasFiniteCoordinates(const Point_3& p)
{
	return std::isfinite(CGAL::to_double(p.x()))
		&& std::isfinite(CGAL::to_double(p.y()))
		&& std::isfinite(CGAL::to_double(p.z()));
}

bool IsDegenerateSoupTriangle(const Point_3& a, const Point_3& b, const Point_3& c)
{
	const double ax = CGAL::to_double(a.x());
	const double ay = CGAL::to_double(a.y());
	const double az = CGAL::to_double(a.z());
	const double bx = CGAL::to_double(b.x());
	const double by = CGAL::to_double(b.y());
	const double bz = CGAL::to_double(b.z());
	const double cx = CGAL::to_double(c.x());
	const double cy = CGAL::to_double(c.y());
	const double cz = CGAL::to_double(c.z());

	const double abx = bx - ax;
	const double aby = by - ay;
	const double abz = bz - az;
	const double acx = cx - ax;
	const double acy = cy - ay;
	const double acz = cz - az;

	const double cross_x = aby * acz - abz * acy;
	const double cross_y = abz * acx - abx * acz;
	const double cross_z = abx * acy - aby * acx;
	const double area_twice_squared = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
	return area_twice_squared <= 1e-24;
}

template <typename PolygonRange>
PolygonSoupStats ComputePolygonSoupStats(const std::vector<Point_3>& points, const PolygonRange& polygons)
{
	PolygonSoupStats stats;
	stats.point_count = points.size();
	stats.triangle_count = polygons.size();
	stats.is_polygon_mesh = PMP::is_polygon_soup_a_polygon_mesh(polygons);

	for (const auto& point : points) {
		if (!HasFiniteCoordinates(point)) {
			++stats.invalid_point_count;
		}
	}

	std::set<TriangleCoordinateKey> triangle_keys;
	for (std::size_t polygon_index = 0; polygon_index < polygons.size(); ++polygon_index) {
		const auto& polygon = polygons[polygon_index];
		if (polygon.size() != 3) {
			++stats.non_triangle_polygon_count;
			continue;
		}

		if (polygon[0] >= points.size() || polygon[1] >= points.size() || polygon[2] >= points.size()) {
			++stats.invalid_index_polygon_count;
			continue;
		}

		const Point_3& p0 = points[polygon[0]];
		const Point_3& p1 = points[polygon[1]];
		const Point_3& p2 = points[polygon[2]];
		if (!HasFiniteCoordinates(p0) || !HasFiniteCoordinates(p1) || !HasFiniteCoordinates(p2)
			|| MakeQuantizedCoordinateKeyForSoup(p0) == MakeQuantizedCoordinateKeyForSoup(p1)
			|| MakeQuantizedCoordinateKeyForSoup(p0) == MakeQuantizedCoordinateKeyForSoup(p2)
			|| MakeQuantizedCoordinateKeyForSoup(p1) == MakeQuantizedCoordinateKeyForSoup(p2)
			|| IsDegenerateSoupTriangle(p0, p1, p2)) {
			++stats.degenerate_triangle_count;
			std::cerr << std::setprecision(17)
				<< "[PrepareOuterBeamSearchNode] Degenerate STL triangle: face=" << polygon_index
				<< ", vertex_indices=(" << polygon[0] << ", " << polygon[1] << ", " << polygon[2] << ")"
				<< ", p0=(" << p0.x() << ", " << p0.y() << ", " << p0.z() << ")"
				<< ", p1=(" << p1.x() << ", " << p1.y() << ", " << p1.z() << ")"
				<< ", p2=(" << p2.x() << ", " << p2.y() << ", " << p2.z() << ")"
				<< std::endl;
		}

		TriangleCoordinateKey key = {
			MakeQuantizedCoordinateKeyForSoup(p0),
			MakeQuantizedCoordinateKeyForSoup(p1),
			MakeQuantizedCoordinateKeyForSoup(p2)
		};
		std::sort(key.begin(), key.end());
		if (!triangle_keys.insert(key).second) {
			++stats.duplicate_triangle_count;
		}
	}

	return stats;
}

void PrintPolygonSoupStats(const std::string& current_file_name, const PolygonSoupStats& stats)
{
	std::cerr << "[PrepareOuterBeamSearchNode] STL polygon soup diagnostics: "
		<< current_file_name
		<< ", points=" << stats.point_count
		<< ", triangles=" << stats.triangle_count
		<< ", non_triangle_polygons=" << stats.non_triangle_polygon_count
		<< ", invalid_index_polygons=" << stats.invalid_index_polygon_count
		<< ", invalid_points=" << stats.invalid_point_count
		<< ", degenerate_triangles=" << stats.degenerate_triangle_count
		<< ", duplicate_triangles=" << stats.duplicate_triangle_count
		<< ", is_polygon_soup_a_polygon_mesh=" << stats.is_polygon_mesh
		<< std::endl;
}

std::string FileStem(const std::string& current_file_name)
{
	const std::size_t slash_pos = current_file_name.find_last_of("/\\");
	std::string name = (slash_pos == std::string::npos) ? current_file_name : current_file_name.substr(slash_pos + 1);
	const std::size_t dot_pos = name.find_last_of('.');
	if (dot_pos != std::string::npos) {
		name = name.substr(0, dot_pos);
	}
	return name;
}

bool CheckCurrentNodeMeshClosedManifold(const std::string& current_file_name, const SurfaceMesh& current_node_mesh)
{
	const bool is_closed = CGAL::is_closed(current_node_mesh);
	std::size_t non_manifold_vertex_count = 0;
	for (auto vertex : current_node_mesh.vertices()) {
		if (PMP::is_non_manifold_vertex(vertex, current_node_mesh)) {
			++non_manifold_vertex_count;
		}
	}

	const bool is_manifold = (non_manifold_vertex_count == 0);
	const bool self_intersects = PMP::does_self_intersect(current_node_mesh);
	const bool bounds_volume = is_closed && !self_intersects && PMP::does_bound_a_volume(current_node_mesh);
	if (is_closed && is_manifold && !self_intersects && bounds_volume) {
		return true;
	}

	std::cerr << "[PrepareOuterBeamSearchNode] current_node_mesh failed validation: "
		<< current_file_name
		<< ", closed=" << is_closed
		<< ", manifold=" << is_manifold
		<< ", non_manifold_vertices=" << non_manifold_vertex_count
		<< ", self_intersects=" << self_intersects
		<< ", bounds_volume=" << bounds_volume
		<< std::endl;

	const std::string file_stem = FileStem(current_file_name);
	std::ofstream boundary_report("model\\" + file_stem + "_mesh_boundary_report.txt");
	if (!boundary_report.is_open()) {
		std::cerr << "[PrepareOuterBeamSearchNode] failed to create boundary diagnostics report for "
			<< current_file_name << std::endl;
		return false;
	}

	boundary_report << std::setprecision(17);
	boundary_report << "source_file: " << current_file_name << '\n';
	boundary_report << "is_closed: " << is_closed << '\n';
	boundary_report << "is_manifold: " << is_manifold << '\n';
	boundary_report << "non_manifold_vertex_count: " << non_manifold_vertex_count << '\n';
	boundary_report << "self_intersects: " << self_intersects << '\n';
	boundary_report << "bounds_volume: " << bounds_volume << '\n';

	std::set<SurfaceMesh::Face_index> boundary_faces;
	if (!is_closed) {
		std::size_t border_edge_count = 0;
		for (auto edge : current_node_mesh.edges()) {
			if (!CGAL::is_border(edge, current_node_mesh)) {
				continue;
			}

			auto h = halfedge(edge, current_node_mesh);
			auto interior_h = CGAL::is_border(h, current_node_mesh) ? opposite(h, current_node_mesh) : h;
			auto s = source(h, current_node_mesh);
			auto t = target(h, current_node_mesh);
			const Point_3 ps = current_node_mesh.point(s);
			const Point_3 pt = current_node_mesh.point(t);
			const auto adjacent_face = face(interior_h, current_node_mesh);
			if (adjacent_face != SurfaceMesh::null_face()) {
				boundary_faces.insert(adjacent_face);
			}

			boundary_report << "border_edge[" << border_edge_count << "]: "
				<< "(" << CGAL::to_double(ps.x()) << ", " << CGAL::to_double(ps.y()) << ", " << CGAL::to_double(ps.z()) << ")"
				<< " -> "
				<< "(" << CGAL::to_double(pt.x()) << ", " << CGAL::to_double(pt.y()) << ", " << CGAL::to_double(pt.z()) << ")"
				<< '\n';
			++border_edge_count;
		}
		boundary_report << "border_edge_count: " << border_edge_count << '\n';
	}
	boundary_report << "boundary_adjacent_face_count: " << boundary_faces.size() << '\n';

	if (!is_manifold) {
		std::size_t non_manifold_index = 0;
		for (auto vertex : current_node_mesh.vertices()) {
			if (!PMP::is_non_manifold_vertex(vertex, current_node_mesh)) {
				continue;
			}

			const Point_3 p = current_node_mesh.point(vertex);
			boundary_report << "non_manifold_vertex[" << non_manifold_index << "]: "
				<< "(" << CGAL::to_double(p.x()) << ", " << CGAL::to_double(p.y()) << ", " << CGAL::to_double(p.z()) << ")"
				<< '\n';
			++non_manifold_index;
		}
		boundary_report << "non_manifold_vertex_count_verified: " << non_manifold_index << '\n';
	}

	if (!boundary_faces.empty()) {
		std::vector<Triangle> boundary_face_triangles;
		boundary_face_triangles.reserve(boundary_faces.size());
		for (auto face_index : boundary_faces) {
			std::vector<Point_3> face_points;
			for (auto vertex_index : CGAL::vertices_around_face(current_node_mesh.halfedge(face_index), current_node_mesh)) {
				face_points.push_back(current_node_mesh.point(vertex_index));
			}
			if (face_points.size() < 3) {
				boundary_report << "boundary_face_with_less_than_3_vertices: face=" << face_index << " vertex_count=" << face_points.size() << '\n';
				continue;
			}

			const Point_3& p0 = face_points[0];
			for (std::size_t i = 1; i + 1 < face_points.size(); ++i) {
				Triangle tri;
				tri._vertices[0] = Vertex(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()), CGAL::to_double(p0.z()));
				tri._vertices[1] = Vertex(CGAL::to_double(face_points[i].x()), CGAL::to_double(face_points[i].y()), CGAL::to_double(face_points[i].z()));
				tri._vertices[2] = Vertex(CGAL::to_double(face_points[i + 1].x()), CGAL::to_double(face_points[i + 1].y()), CGAL::to_double(face_points[i + 1].z()));
				tri.normal = Vertex(0, 0, 0);
				boundary_face_triangles.push_back(tri);
			}
		}
		if (!boundary_face_triangles.empty()) {
			Katana::Instance().stl.saveStl("model\\" + file_stem + "_boundary_faces.stl", boundary_face_triangles);
			std::cerr << "[PrepareOuterBeamSearchNode] boundary-adjacent faces saved to model\\"
				<< file_stem << "_boundary_faces.stl" << std::endl;
		}
	}

	boundary_report.close();
	std::cerr << "[PrepareOuterBeamSearchNode] boundary diagnostics saved to model\\"
		<< file_stem << "_mesh_boundary_report.txt" << std::endl;
	return false;
}

bool LoadCurrentNodeMeshFromPolygonSoupForDiagnostics(const std::string& current_file_name, SurfaceMesh& current_node_mesh)
{
	std::vector<Point_3> points;
	std::vector<std::vector<std::size_t>> polygons;
	if (!CGAL::IO::read_STL(current_file_name, points, polygons, CGAL::parameters::verbose(true))) {
		std::cerr << "[PrepareOuterBeamSearchNode] Failed to read STL polygon soup: "
			<< current_file_name << std::endl;
		return false;
	}

	const PolygonSoupStats stats = ComputePolygonSoupStats(points, polygons);
	PrintPolygonSoupStats(current_file_name, stats);

	if (stats.point_count == 0 || stats.triangle_count == 0) {
		std::cerr << "[PrepareOuterBeamSearchNode] Skip polygon_soup_to_polygon_mesh because the STL soup is empty: "
			<< current_file_name << std::endl;
		return false;
	}

	if (!stats.is_polygon_mesh) {
		std::cerr << "[PrepareOuterBeamSearchNode] Skip polygon_soup_to_polygon_mesh because the soup is not a valid polygon mesh: "
			<< current_file_name << std::endl;
		return false;
	}

	current_node_mesh.clear();
	PMP::polygon_soup_to_polygon_mesh(points, polygons, current_node_mesh);
	std::cerr << "[PrepareOuterBeamSearchNode] polygon_soup_to_polygon_mesh result: "
		<< current_file_name
		<< ", vertices=" << current_node_mesh.number_of_vertices()
		<< ", faces=" << current_node_mesh.number_of_faces()
		<< std::endl;

	return true;
}

} // namespace

bool LoadAndCheckCurrentNodeMesh(const std::string& current_file_name, SurfaceMesh& current_node_mesh)
{
	if (!CGAL::IO::read_polygon_mesh(current_file_name, current_node_mesh)) {
		std::cerr << "[PrepareOuterBeamSearchNode] Failed to read polygon mesh: "
			<< current_file_name << std::endl;

		if (!LoadCurrentNodeMeshFromPolygonSoupForDiagnostics(current_file_name, current_node_mesh)) {
			return false;
		}

		std::cerr << "[PrepareOuterBeamSearchNode] Continue with mesh converted from polygon soup: "
			<< current_file_name
			<< std::endl;
	}

	return CheckCurrentNodeMeshClosedManifold(current_file_name, current_node_mesh);
}

} // namespace vasco::mesh_validation
