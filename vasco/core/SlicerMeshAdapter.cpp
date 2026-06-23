#include "vasco/core/SlicerMeshAdapter.h"
#include "vasco/core/Types.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <map>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/helpers.h>

namespace vasco::slicer_mesh_adapter {
namespace {

namespace PMP = CGAL::Polygon_mesh_processing;
using Tri3 = vasco::core::Tri3;

double SquaredDistanceForSlicerPoint(const vasco::Slicer::Vec3& a, const vasco::Slicer::Vec3& b)
{
	const double dx = a[0] - b[0];
	const double dy = a[1] - b[1];
	const double dz = a[2] - b[2];
	return dx * dx + dy * dy + dz * dz;
}

bool IsDegenerateSlicerTriangle(const vasco::Slicer& slicer, const Tri3& tri, double eps = 1e-12)
{
	if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0
		|| tri[0] >= static_cast<int>(slicer.positions.size())
		|| tri[1] >= static_cast<int>(slicer.positions.size())
		|| tri[2] >= static_cast<int>(slicer.positions.size())) {
		return true;
	}

	const auto& a = slicer.positions[tri[0]];
	const auto& b = slicer.positions[tri[1]];
	const auto& c = slicer.positions[tri[2]];
	if (SquaredDistanceForSlicerPoint(a, b) <= eps
		|| SquaredDistanceForSlicerPoint(a, c) <= eps
		|| SquaredDistanceForSlicerPoint(b, c) <= eps) {
		return true;
	}

	const double abx = b[0] - a[0];
	const double aby = b[1] - a[1];
	const double abz = b[2] - a[2];
	const double acx = c[0] - a[0];
	const double acy = c[1] - a[1];
	const double acz = c[2] - a[2];
	const double cross_x = aby * acz - abz * acy;
	const double cross_y = abz * acx - abx * acz;
	const double cross_z = abx * acy - aby * acx;
	return cross_x * cross_x + cross_y * cross_y + cross_z * cross_z <= eps;
}

bool IsTriangleOnCutPlane(const vasco::Slicer& slicer, const Tri3& tri, const std::vector<double>& cut_plane_z_values, double eps = 1e-4)
{
	if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0
		|| tri[0] >= static_cast<int>(slicer.positions.size())
		|| tri[1] >= static_cast<int>(slicer.positions.size())
		|| tri[2] >= static_cast<int>(slicer.positions.size())) {
		return false;
	}

	for (double z : cut_plane_z_values) {
		if (std::abs(slicer.positions[tri[0]][2] - z) <= eps
			&& std::abs(slicer.positions[tri[1]][2] - z) <= eps
			&& std::abs(slicer.positions[tri[2]][2] - z) <= eps) {
			return true;
		}
	}
	return false;
}

void CollectCutPatchFaces(
	const vasco::Slicer& slicer,
	const std::vector<int>& id_contact_faces,
	std::vector<bool>& selected_faces)
{
	selected_faces.assign(slicer.triangles.size(), false);
	std::vector<std::vector<int>> vertex_to_faces(slicer.positions.size());
	for (int face_id = 0; face_id < static_cast<int>(slicer.triangles.size()); ++face_id) {
		const auto& tri = slicer.triangles[face_id];
		for (int k = 0; k < 3; ++k) {
			if (tri[k] >= 0 && tri[k] < static_cast<int>(vertex_to_faces.size())) {
				vertex_to_faces[tri[k]].push_back(face_id);
			}
		}
	}

	std::vector<int> frontier;
	for (int face_id : id_contact_faces) {
		if (face_id >= 0 && face_id < static_cast<int>(selected_faces.size())) {
			selected_faces[face_id] = true;
			frontier.push_back(face_id);
		}
	}

	for (int face_id : frontier) {
		const auto& tri = slicer.triangles[face_id];
		for (int k = 0; k < 3; ++k) {
			if (tri[k] < 0 || tri[k] >= static_cast<int>(vertex_to_faces.size())) {
				continue;
			}
			for (int neighbor_face_id : vertex_to_faces[tri[k]]) {
				selected_faces[neighbor_face_id] = true;
			}
		}
	}
}

bool ConvertSlicerToSurfaceMesh(
	const vasco::Slicer& slicer,
	SurfaceMesh& mesh,
	std::vector<SurfaceMesh::Face_index>& old_face_to_new_face)
{
	mesh.clear();
	old_face_to_new_face.assign(slicer.triangles.size(), SurfaceMesh::null_face());
	std::vector<SurfaceMesh::Vertex_index> vertex_map;
	vertex_map.reserve(slicer.positions.size());
	for (const auto& p : slicer.positions) {
		vertex_map.push_back(mesh.add_vertex(Point_3(p[0], p[1], p[2])));
	}

	std::size_t skipped_degenerate_faces = 0;
	std::size_t add_face_failed_count = 0;
	for (int face_id = 0; face_id < static_cast<int>(slicer.triangles.size()); ++face_id) {
		const auto& tri = slicer.triangles[face_id];
		if (IsDegenerateSlicerTriangle(slicer, tri)) {
			++skipped_degenerate_faces;
			continue;
		}
		auto face = mesh.add_face(vertex_map[tri[0]], vertex_map[tri[1]], vertex_map[tri[2]]);
		if (face == SurfaceMesh::null_face()) {
			++add_face_failed_count;
			continue;
		}
		old_face_to_new_face[face_id] = face;
	}

	if (skipped_degenerate_faces != 0) {
		std::cerr << "[HybridManufacturing::CutMesh] Skip degenerate faces before local remesh: "
			<< skipped_degenerate_faces << std::endl;
	}
	if (add_face_failed_count != 0) {
		std::cerr << "[HybridManufacturing::CutMesh] Skip local remesh because SurfaceMesh rejected faces: "
			<< add_face_failed_count << std::endl;
		return false;
	}

	return mesh.number_of_faces() != 0;
}

void ConvertSurfaceMeshToSlicer(const SurfaceMesh& mesh, vasco::Slicer& slicer)
{
	slicer.clear();
	std::map<SurfaceMesh::Vertex_index, int> vertex_index_map;
	for (auto vertex : mesh.vertices()) {
		const auto& point = mesh.point(vertex);
		const int index = static_cast<int>(slicer.positions.size());
		slicer.positions.push_back({ CGAL::to_double(point.x()), CGAL::to_double(point.y()), CGAL::to_double(point.z()) });
		vertex_index_map[vertex] = index;
	}

	for (auto face : mesh.faces()) {
		Tri3 tri{};
		int vertex_count = 0;
		for (auto halfedge : CGAL::halfedges_around_face(mesh.halfedge(face), mesh)) {
			if (vertex_count >= 3) {
				break;
			}
			tri[vertex_count++] = vertex_index_map[mesh.target(halfedge)];
		}
		if (vertex_count == 3 && !IsDegenerateSlicerTriangle(slicer, tri)) {
			slicer.triangles.push_back(tri);
		}
	}
}

} // namespace

void RemeshCutPatchBeforeSaving(
	vasco::Slicer& all_slicer,
	std::vector<int>& id_contact_faces,
	const std::vector<double>& cut_plane_z_values,
	double target_edge_length)
{
	if (id_contact_faces.empty() || all_slicer.triangles.empty()) {
		return;
	}

	std::vector<bool> selected_old_faces;
	CollectCutPatchFaces(all_slicer, id_contact_faces, selected_old_faces);

	SurfaceMesh mesh;
	std::vector<SurfaceMesh::Face_index> old_face_to_new_face;
	if (!ConvertSlicerToSurfaceMesh(all_slicer, mesh, old_face_to_new_face)) {
		return;
	}

	std::vector<SurfaceMesh::Face_index> selected_faces;
	for (int face_id = 0; face_id < static_cast<int>(selected_old_faces.size()); ++face_id) {
		if (selected_old_faces[face_id] && old_face_to_new_face[face_id] != SurfaceMesh::null_face()) {
			selected_faces.push_back(old_face_to_new_face[face_id]);
		}
	}
	if (selected_faces.empty()) {
		return;
	}

	try {
		PMP::isotropic_remeshing(
			selected_faces,
			target_edge_length,
			mesh,
			CGAL::parameters::number_of_iterations(2)
				.protect_constraints(true));
	}
	catch (const std::exception& e) {
		std::cerr << "[HybridManufacturing::CutMesh] Local cut patch remesh failed: "
			<< e.what() << std::endl;
		return;
	}

	ConvertSurfaceMeshToSlicer(mesh, all_slicer);

	id_contact_faces.clear();
	for (int face_id = 0; face_id < static_cast<int>(all_slicer.triangles.size()); ++face_id) {
		if (IsTriangleOnCutPlane(all_slicer, all_slicer.triangles[face_id], cut_plane_z_values)) {
			id_contact_faces.push_back(face_id);
		}
	}

	std::cerr << "[HybridManufacturing::CutMesh] Local cut patch remesh finished: target_edge_length="
		<< target_edge_length
		<< ", contact_faces=" << id_contact_faces.size()
		<< ", vertices=" << all_slicer.positions.size()
		<< ", faces=" << all_slicer.triangles.size()
		<< std::endl;
}

} // namespace vasco::slicer_mesh_adapter
