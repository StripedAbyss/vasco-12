#include "HybridManufacturing.h"
#include<stdlib.h>
#include<igl/readOBJ.h>
#include<igl/writeOBJ.h>
#include <array>
#include <exception>
#include <iomanip>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/IO/STL.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

namespace
{
	using ContactTriangulationTraits = CGAL::Exact_predicates_inexact_constructions_kernel;
	using ContactPoint2 = ContactTriangulationTraits::Point_2;
	using ContactVb = CGAL::Triangulation_vertex_base_2<ContactTriangulationTraits>;
	template <typename GT, typename Fb = CGAL::Constrained_triangulation_face_base_2<GT>>
	class ContactFaceBase : public Fb {
	public:
		bool is_in_domain() const { return in_domain_; }
		void set_in_domain(bool v) { in_domain_ = v; }

		template <typename TDS2>
		struct Rebind_TDS {
			using Fb2 = typename Fb::template Rebind_TDS<TDS2>::Other;
			using Other = ContactFaceBase<GT, Fb2>;
		};

		ContactFaceBase() : Fb() {}
		ContactFaceBase(typename Fb::Vertex_handle v0, typename Fb::Vertex_handle v1, typename Fb::Vertex_handle v2) : Fb(v0, v1, v2) {}
		ContactFaceBase(typename Fb::Vertex_handle v0, typename Fb::Vertex_handle v1, typename Fb::Vertex_handle v2,
			typename Fb::Face_handle n0, typename Fb::Face_handle n1, typename Fb::Face_handle n2)
			: Fb(v0, v1, v2, n0, n1, n2) {}

	private:
		bool in_domain_ = false;
	};

	using ContactFb = ContactFaceBase<ContactTriangulationTraits>;
	using ContactTds = CGAL::Triangulation_data_structure_2<ContactVb, ContactFb>;
	using ContactCdt = CGAL::Constrained_Delaunay_triangulation_2<ContactTriangulationTraits, ContactTds>;
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
		if (!is_closed || !is_manifold) {
			std::cerr << "[PrepareOuterBeamSearchNode] current_node_mesh is not a closed manifold mesh: "
				<< current_file_name
				<< ", closed=" << is_closed
				<< ", manifold=" << is_manifold
				<< ", non_manifold_vertices=" << non_manifold_vertex_count
				<< std::endl;
			const std::string file_stem = [&current_file_name]() {
				const std::size_t slash_pos = current_file_name.find_last_of("/\\");
				std::string name = (slash_pos == std::string::npos) ? current_file_name : current_file_name.substr(slash_pos + 1);
				const std::size_t dot_pos = name.find_last_of('.');
				if (dot_pos != std::string::npos) {
					name = name.substr(0, dot_pos);
				}
				return name;
			}();

			std::ofstream boundary_report("model\\" + file_stem + "_mesh_boundary_report.txt");
			if (boundary_report.is_open()) {
				boundary_report << std::setprecision(17);
				boundary_report << "source_file: " << current_file_name << '\n';
				boundary_report << "is_closed: " << is_closed << '\n';
				boundary_report << "is_manifold: " << is_manifold << '\n';
				boundary_report << "non_manifold_vertex_count: " << non_manifold_vertex_count << '\n';

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
			}
			else {
				std::cerr << "[PrepareOuterBeamSearchNode] failed to create boundary diagnostics report for "
					<< current_file_name << std::endl;
			}
			return false;
		}

		return true;
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

	double SquaredDistanceForSlicerPoint(const Slicer_2::Vec3& a, const Slicer_2::Vec3& b)
	{
		const double dx = a[0] - b[0];
		const double dy = a[1] - b[1];
		const double dz = a[2] - b[2];
		return dx * dx + dy * dy + dz * dz;
	}

	bool IsDegenerateSlicerTriangle(const Slicer_2& slicer, const TRiangle& tri, double eps = 1e-12)
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

	bool IsTriangleOnCutPlane(const Slicer_2& slicer, const TRiangle& tri, const std::vector<double>& cut_plane_z_values, double eps = 1e-4)
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
		const Slicer_2& slicer,
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
		const Slicer_2& slicer,
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

	void ConvertSurfaceMeshToSlicer(const SurfaceMesh& mesh, Slicer_2& slicer)
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
			TRiangle tri{};
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

	void RemeshCutPatchBeforeSaving(
		Slicer_2& all_slicer,
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

	template <typename Point2>
	void AppendUniquePoint(std::vector<Point2>& points, const Point2& point)
	{
		if (points.empty()) {
			points.push_back(point);
			return;
		}
		const double eps = 1e-12;
		const auto& last = points.back();
		if (std::abs(CGAL::to_double(last.x() - point.x())) < eps && std::abs(CGAL::to_double(last.y() - point.y())) < eps) {
			return;
		}
		points.push_back(point);
	}

	inline double SignedArea2D(const Point_3& a, const Point_3& b, const Point_3& c)
	{
		return (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x());
	}

	struct QuantizedPointKey {
		long long x = 0;
		long long y = 0;
		long long z = 0;
		bool operator<(const QuantizedPointKey& other) const
		{
			if (x != other.x) return x < other.x;
			if (y != other.y) return y < other.y;
			return z < other.z;
		}
	};

	inline QuantizedPointKey MakeQuantizedKey(const Point_3& p, double scale = 1e9)
	{
		return {
			static_cast<long long>(std::llround(CGAL::to_double(p.x()) * scale)),
			static_cast<long long>(std::llround(CGAL::to_double(p.y()) * scale)),
			static_cast<long long>(std::llround(CGAL::to_double(p.z()) * scale))
		};
	}

	inline int GetOrCreateSlicerPositionIndex(
		std::map<QuantizedPointKey, int>& point_index_map,
		vasco::Slicer& slicer,
		const Point_3& point)
	{
		auto key = MakeQuantizedKey(point);
		auto it = point_index_map.find(key);
		if (it != point_index_map.end()) {
			return it->second;
		}
		int index = static_cast<int>(slicer.positions.size());
		slicer.positions.push_back({ CGAL::to_double(point.x()), CGAL::to_double(point.y()), CGAL::to_double(point.z()) });
		point_index_map.emplace(key, index);
		return index;
	}

	void AddConstraintLoop(ContactCdt& cdt, const std::vector<int>& ring, const Slicer_2& slicer)
	{
		if (ring.size() < 3) return;
		std::vector<ContactCdt::Vertex_handle> handles;
		handles.reserve(ring.size());
		for (int idx : ring) {
			const auto& pos = slicer.positions[idx];
			handles.push_back(cdt.insert(ContactPoint2(pos[0], pos[1])));
		}
		for (size_t i = 0; i < handles.size(); ++i) {
			cdt.insert_constraint(handles[i], handles[(i + 1) % handles.size()]);
		}
	}

	double ComputeAverageBoundaryEdgeLength(const std::vector<int>& ring, const Slicer_2& slicer)
	{
		if (ring.size() < 2) {
			return 1.0;
		}
		double sum = 0.0;
		for (size_t i = 0; i < ring.size(); ++i) {
			const auto& a = slicer.positions[ring[i]];
			const auto& b = slicer.positions[ring[(i + 1) % ring.size()]];
			double dx = a[0] - b[0];
			double dy = a[1] - b[1];
			sum += std::sqrt(dx * dx + dy * dy);
		}
		return std::max(sum / static_cast<double>(ring.size()), 1e-6);
	}

	std::vector<TRiangle> TriangulateContactFacesWithCDT(
		const std::vector<int>& outer_ring,
		const std::vector<std::vector<int>>& hole_rings,
		Slicer_2& slicer,
		std::vector<int>& id_contact_faces,
		std::map<QuantizedPointKey, int>& point_index_map)
	{
		std::vector<TRiangle> generated;
		if (outer_ring.size() < 3) {
			return generated;
		}

		ContactCdt cdt;
		AddConstraintLoop(cdt, outer_ring, slicer);
		for (const auto& hole : hole_rings) {
			AddConstraintLoop(cdt, hole, slicer);
		}

		CGAL::mark_domain_in_triangulation(cdt);

		double target_edge = ComputeAverageBoundaryEdgeLength(outer_ring, slicer);
		typedef CGAL::Delaunay_mesh_size_criteria_2<ContactCdt> Criteria;
		Criteria criteria(0.125, target_edge);
		CGAL::refine_Delaunay_mesh_2(cdt, criteria, true);

		std::map<ContactCdt::Vertex_handle, int> vertex_index_map;
		for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
			const auto p = vit->point();
			Point_3 point(CGAL::to_double(p.x()), CGAL::to_double(p.y()), slicer.positions[outer_ring.front()][2]);
			vertex_index_map[vit] = GetOrCreateSlicerPositionIndex(point_index_map, slicer, point);
		}

		for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
			if (!fit->is_in_domain()) {
				continue;
			}
			TRiangle tri{};
			tri[0] = vertex_index_map[fit->vertex(0)];
			tri[1] = vertex_index_map[fit->vertex(1)];
			tri[2] = vertex_index_map[fit->vertex(2)];
			const auto& p0 = slicer.positions[tri[0]];
			const auto& p1 = slicer.positions[tri[1]];
			const auto& p2 = slicer.positions[tri[2]];
			Point_3 a(p0[0], p0[1], p0[2]);
			Point_3 b(p1[0], p1[1], p1[2]);
			Point_3 c(p2[0], p2[1], p2[2]);
			if (SignedArea2D(a, b, c) < 0.0) {
				std::swap(tri[1], tri[2]);
			}
			generated.push_back(tri);
			slicer.triangles.push_back(tri);
			id_contact_faces.push_back(static_cast<int>(slicer.triangles.size() - 1));
		}
		return generated;
	}
}

namespace
{
	inline void PrepareToolForCollision(cutter& tool)
	{
		tool.cylinder_height_threshold = tool.cylinder_height + tool.ball_r;
		tool.carriage_check_radius_sq = (tool.carriage_r + 5.0) * (tool.carriage_r + 5.0);
		tool.cylinder_check_radius_sq = (tool.cylinder_r + 5.0) * (tool.cylinder_r + 5.0);
		tool.cylinder_r_sq = tool.cylinder_r * tool.cylinder_r;
		tool.carriage_r_sq = tool.carriage_r * tool.carriage_r;
		tool.total_height = tool.cylinder_height + tool.ball_r + tool.carriage_height;
	}

	inline Eigen::Vector3d ToVector3(const Eigen::MatrixXd& vec)
	{
		return { vec(0, 0), vec(1, 0), vec(2, 0) };
	}

	inline Eigen::Vector3d ComputeToolCenter(const Eigen::Vector3d& point, const Eigen::Vector3d& normal, double radius)
	{
		return point + radius * normal;
	}

	inline Vertex LocalComputeIntersection(const Vertex& a, const Vertex& b, float z)
	{
		float t = (z - a.z) / (b.z - a.z);
		Vertex intersection = { a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), z };
		return intersection;
	}

	inline Vertex LocalCalTriNormal(const Vertex& ver1, const Vertex& ver2, const Vertex& ver3)
	{
		float temp1[3], temp2[3], normal[3];
		temp1[0] = ver2.x - ver1.x;
		temp1[1] = ver2.y - ver1.y;
		temp1[2] = ver2.z - ver1.z;
		temp2[0] = ver3.x - ver2.x;
		temp2[1] = ver3.y - ver2.y;
		temp2[2] = ver3.z - ver2.z;

		normal[0] = temp1[1] * temp2[2] - temp1[2] * temp2[1];
		normal[1] = -(temp1[0] * temp2[2] - temp1[2] * temp2[0]);
		normal[2] = temp1[0] * temp2[1] - temp1[1] * temp2[0];

		float length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
		if (length == 0.0f) {
			length = 1.0f;
		}
		normal[0] /= length;
		normal[1] /= length;
		normal[2] /= length;
		return Vertex(normal[0], normal[1], normal[2]);
	}

	inline Segment LocalComputeSegment(Triangle& t, float z)
	{
		Vertex** vs = t.vertices;
		Segment segment;
		segment.neighbours[0] = NULL;
		segment.neighbours[1] = NULL;
		segment.orderIndex = -1;

		assert(vs[0]->z <= vs[1]->z);
		assert(vs[1]->z <= vs[2]->z);
		assert(z >= vs[0]->z);
		assert(z <= vs[2]->z);

		if (z < vs[1]->z) {
			segment.vertices[0] = LocalComputeIntersection(*vs[0], *vs[1], z);
			segment.vertices[1] = LocalComputeIntersection(*vs[0], *vs[2], z);
		}
		else {
			segment.vertices[0] = LocalComputeIntersection(*vs[1], *vs[2], z);
			segment.vertices[1] = LocalComputeIntersection(*vs[0], *vs[2], z);
		}

		segment.triangle_points.push_back(*t.vertices[0]);
		segment.triangle_points.push_back(*t.vertices[1]);
		segment.triangle_points.push_back(*t.vertices[2]);

		t.normal = LocalCalTriNormal(*vs[0], *vs[1], *vs[2]);
		assert(t.normal.z < 0.995 && t.normal.z > -0.995);
		Vertex n = t.normal;
		n.z = 0;
		n = n.normalize();
		segment.normal = n;
		return segment;
	}

	inline void LocalUnifySegmentVertices(std::vector<Segment>& segments, std::map<Vertex, std::vector<Segment*>>& segmentsByVertex)
	{
		for (unsigned int i = 0; i < segments.size(); i++) {
			Segment& s = segments[i];
			for (int j = 0; j < 2; j++) {
				Vertex& v = s.vertices[j];
				segmentsByVertex[v].push_back(&s);
			}
		}
	}

	inline bool LocalCheckInside(Point pt, Point* pgn_begin, Point* pgn_end)
	{
		switch (CGAL::bounded_side_2(pgn_begin, pgn_end, pt, K())) {
		case CGAL::ON_BOUNDED_SIDE:
			return true;
		case CGAL::ON_BOUNDARY:
		case CGAL::ON_UNBOUNDED_SIDE:
		default:
			return false;
		}
	}
}

void HybridManufacturing::BuildSurfaceMeshSlices(std::vector<SurfaceMeshSliceData>& slices) const
{
	slices.clear();
	const double layer_height_value = Katana::Instance().config.get("layer_height");
	if (layer_height_value <= 0) {
		return;
	}

	//std::vector<Eigen::Vector3d> face_normals;
	//face_normals.reserve(current_node_mesh_rotated.number_of_faces());
	//int face_id = 0;
	//for (auto face : current_node_mesh_rotated.faces()) {
	//	auto halfedge = current_node_mesh_rotated.halfedge(face);
	//	auto v0 = current_node_mesh_rotated.target(halfedge);
	//	auto v1 = current_node_mesh_rotated.target(current_node_mesh_rotated.next(halfedge));
	//	auto v2 = current_node_mesh_rotated.target(current_node_mesh_rotated.next(current_node_mesh_rotated.next(halfedge)));
	//	const auto& p0 = current_node_mesh_rotated.point(v0);
	//	const auto& p1 = current_node_mesh_rotated.point(v1);
	//	const auto& p2 = current_node_mesh_rotated.point(v2);
	//	Eigen::Vector3d e1(p1.x() - p0.x(), p1.y() - p0.y(), p1.z() - p0.z());
	//	Eigen::Vector3d e2(p2.x() - p0.x(), p2.y() - p0.y(), p2.z() - p0.z());
	//	Eigen::Vector3d normal = e1.cross(e2);
	//	if (normal.norm() > 0.0) {
	//		normal.normalize();
	//	}
	//	face_normals.push_back(normal);
	//}
	//slice_data.face_normals = face_normals;


	std::vector<SurfaceMesh::Vertex_index> vi_z_sorted;
	vi_z_sorted.reserve(current_node_mesh_rotated.number_of_vertices());
	for (auto v : current_node_mesh_rotated.vertices()) {
		vi_z_sorted.push_back(v);
	}
	std::sort(vi_z_sorted.begin(), vi_z_sorted.end(), [&](const SurfaceMesh::Vertex_index& a, const SurfaceMesh::Vertex_index& b) {
		return current_node_mesh_rotated.point(a).z() < current_node_mesh_rotated.point(b).z();
		});

	std::map<SurfaceMesh::Face_index, int> activeTriangles;

	double min_z = current_node_mesh_rotated.point(vi_z_sorted.front()).z();
	double layer_z = min_z + layer_height_value;
	int layer_height_adjust_count = 0;
	for (auto vi : vi_z_sorted) {
		double vi_z = current_node_mesh_rotated.point(vi).z();
		const auto& point = current_node_mesh_rotated.point(vi);
		//std::cout << "vi: " << vi << " z: " << vi_z << std::endl;

		while (layer_z <= vi_z) {
			if (vi_z - layer_z < layer_vertex_threshold) {
				layer_z += small_layer_height;
				++layer_height_adjust_count;
				std::cout << "Adjusting layer height for vertex at z = " << vi_z << " layer_height_adjust_count: " << layer_height_adjust_count << std::endl;
				continue;
			}


			SurfaceMeshSliceData slice_data;
			if (!BuildSurfaceMeshSingleSlice(activeTriangles, slice_data, layer_z)) {
				layer_z += small_layer_height;
				++layer_height_adjust_count;
				std::cout << "Adjusting layer height for vertex at z = " << vi_z << " layer_height_adjust_count: " << layer_height_adjust_count << std::endl;
				continue;
			}

			// build slice successfully
			layer_height_adjust_count = 0;
			if (slice_data.contour_points.empty()) {
				std::cout << "empty slice at z = " << layer_z << std::endl;
			}
			else {
				slices.push_back(std::move(slice_data));
			}
			layer_z += layer_height_value;
		}

		// Update activeTriangles
		auto hf_vi = current_node_mesh_rotated.halfedge(vi);
		for (auto hf : current_node_mesh_rotated.halfedges_around_target(hf_vi)) {
			auto face_index = current_node_mesh_rotated.face(hf);
			if (face_index == SurfaceMesh::null_face()) {
				std::cout << "[HybridManufacturing::BuildSurfaceMeshSlices] Vertex " << vi << " has a halfedge with null face." << std::endl;
				continue;
			}
			auto active_it = activeTriangles.find(face_index);
			if (active_it == activeTriangles.end()) {
				activeTriangles[face_index] = 1;
			}
			else {
				active_it->second += 1;
				if (active_it->second == 3) {
					activeTriangles.erase(active_it);
				}
			}
		}

	}
	if (activeTriangles.size() > 0) {
		std::cout << "[HybridManufacturing::BuildSurfaceMeshSlices] Warning: " << activeTriangles.size() << " active triangles remain after processing all vertices." << std::endl;
	}
}

bool HybridManufacturing::BuildSurfaceMeshSingleSlice(std::map<SurfaceMesh::Face_index, int>& activeTriangles, SurfaceMeshSliceData& slice_data, double layer_z) const
{
	if (activeTriangles.empty()) {
		return true;
	}

	std::vector<bool> face_visited(activeTriangles.size(), false);
	std::vector<SurfaceMesh::Face_index> active_face_indices;
	active_face_indices.reserve(activeTriangles.size());
	std::unordered_map<SurfaceMesh::Face_index, int> face_index_to_active_index;
	for (auto active_triangle : activeTriangles) {
		const auto& face_index = active_triangle.first;
		active_face_indices.push_back(face_index);
		face_index_to_active_index[face_index] = active_face_indices.size() - 1;
	}

	for (size_t i = 0; i < active_face_indices.size(); ++i) {
		if (face_visited[i]) {
			continue;
		}
		const auto start_face_index = active_face_indices[i];

		// Check if the face has vertices close to the slicing plane
		auto hf = current_node_mesh_rotated.halfedge(start_face_index);
		for (int j = 0; j < 3; j++) {
			auto v = current_node_mesh_rotated.target(hf);
			const auto& p = current_node_mesh_rotated.point(v);
			if (std::abs(p.z() - layer_z) < layer_vertex_threshold) {
				return false;
			}
			hf = current_node_mesh_rotated.next(hf);
		}
		if (hf != current_node_mesh_rotated.halfedge(start_face_index)) {
			std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << start_face_index << " does not have 3 vertices." << std::endl;
			return false;
		}

		// Find the intersection of the face with the slicing plane
		bool intersected = false;
		for (int j = 0; j < 3; j++) {
			auto v0 = current_node_mesh_rotated.source(hf);
			auto v1 = current_node_mesh_rotated.target(hf);
			const auto& p0 = current_node_mesh_rotated.point(v0);
			const auto& p1 = current_node_mesh_rotated.point(v1);
			if (!(p0.z() < layer_z && p1.z() > layer_z)) {
				hf = current_node_mesh_rotated.next(hf);
				continue;
			}
			intersected = true;
			Point_3 intersection_point = p0 + (layer_z - p0.z()) / (p1.z() - p0.z()) * (p1 - p0);
			slice_data.contour_points.push_back(Polyline_type{ intersection_point });
			slice_data.contour_face_ids.push_back({ start_face_index });
			slice_data.contour_contain_next_id.push_back(-1);
			slice_data.layer_z = layer_z;

			break;
		}

		if (!intersected) {
			std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << start_face_index << " does not intersect with the slicing plane. Vertices: ";
			if (start_face_index == SurfaceMesh::null_face()) {
				std::cout << "null face" << std::endl;
				return false;
			}
			auto debug_hf = current_node_mesh_rotated.halfedge(start_face_index);
			for (int k = 0; k < 3; ++k) {
				const auto& p = current_node_mesh_rotated.point(current_node_mesh_rotated.target(debug_hf));
				std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
				if (k < 2) {
					std::cout << ", ";
				}
				debug_hf = current_node_mesh_rotated.next(debug_hf);
			}
			std::cout << std::endl;
			std::cout << "layer_z: " << layer_z << std::endl;
			return false;
		}

		face_visited[i] = true;

		const auto start_halfedge = hf;

		hf = current_node_mesh_rotated.opposite(hf);
		while (current_node_mesh_rotated.face(hf) != start_face_index) {
			// Find the intersection of the face with the slicing plane
			bool intersected = false;
			for (int j = 0; j < 2; j++) {
				hf = current_node_mesh_rotated.next(hf);
				auto v0 = current_node_mesh_rotated.source(hf);
				auto v1 = current_node_mesh_rotated.target(hf);
				const auto& p0 = current_node_mesh_rotated.point(v0);
				const auto& p1 = current_node_mesh_rotated.point(v1);
				if (!(p0.z() < layer_z && p1.z() > layer_z)) {
					continue;
				}
				intersected = true;
				Point_3 intersection_point = p0 + (layer_z - p0.z()) / (p1.z() - p0.z()) * (p1 - p0);
				auto face_index = current_node_mesh_rotated.face(hf);
				auto active_index_it = face_index_to_active_index.find(face_index);

				if (active_index_it == face_index_to_active_index.end()) {
					std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << face_index << " is not in activeTriangles." << std::endl;
					return false;
				}
				if (face_visited[active_index_it->second]) {
					std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << face_index << " has already been visited." << std::endl;
					return false;
				}
				face_visited[active_index_it->second] = true;


				slice_data.contour_points.back().push_back(intersection_point);
				slice_data.contour_face_ids.back().push_back(face_index);

				break;
			}

			if (!intersected) {
				auto face_index = current_node_mesh_rotated.face(hf);
				std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << face_index << " does not intersect with the slicing plane." << std::endl;
				if (face_index == SurfaceMesh::null_face()) {
					std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Halfedge " << hf << " has null face." << std::endl;
					return false;
				}
				auto debug_hf = current_node_mesh_rotated.halfedge(face_index);
				if (debug_hf == SurfaceMesh::null_halfedge()) {
					std::cout << "[HybridManufacturing::BuildSurfaceMeshSingleSlice] Face " << face_index << " has null halfedge." << std::endl;
					return false;
				}
				for (int k = 0; k < 3; ++k) {
					const auto& p = current_node_mesh_rotated.point(current_node_mesh_rotated.target(debug_hf));
					std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
					if (k < 2) {
						std::cout << ", ";
					}
					debug_hf = current_node_mesh_rotated.next(debug_hf);
				}
				std::cout << std::endl;
				std::cout << "layer_z: " << layer_z << std::endl;
				return false;
				return false;
			}

			hf = current_node_mesh_rotated.opposite(hf);
		}
	}

	//TODO: add contain relationship between contours
	return true;
}

namespace
{
	inline double PointDistance2D(const Point_3& a, const Point_3& b)
	{
		double dx = a.x() - b.x();
		double dy = a.y() - b.y();
		return std::sqrt(dx * dx + dy * dy);
	}

	inline bool PointsClose(const Point_3& a, const Point_3& b, double tol)
	{
		return PointDistance2D(a, b) <= tol;
	}
}

HybridManufacturing::HybridManufacturing(std::string file_name, std::string suf, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd N)
{
	this->V = V;
	this->F = F;
	this->Normals = N;
	this->file_name = file_name;
	this->suf = suf;
	InitializeSurfaceMeshFromVF();
	InitializePolyscope();
}

HybridManufacturing::~HybridManufacturing()
{
}

void HybridManufacturing::SetContactFaceTriangulationMode(ContactFaceTriangulationMode mode)
{
	contact_face_triangulation_mode = mode;
}

HybridManufacturing::ContactFaceTriangulationMode HybridManufacturing::GetContactFaceTriangulationMode() const
{
	return contact_face_triangulation_mode;
}

void HybridManufacturing::InitializeSurfaceMeshFromVF()
{
	input_mesh.clear();
	if (V.rows() == 0 || F.rows() == 0) {
		return;
	}

	std::vector<SurfaceMesh::Vertex_index> vertices;
	vertices.reserve(V.rows());
	for (int i = 0; i < V.rows(); ++i) {
		vertices.push_back(input_mesh.add_vertex(Point_3(V(i, 0), V(i, 1), V(i, 2))));
	}

	for (int i = 0; i < F.rows(); ++i) {
		input_mesh.add_face(vertices[F(i, 0)], vertices[F(i, 1)], vertices[F(i, 2)]);
	}
}

std::vector<TRiangle> HybridManufacturing::FilterSurfaceRemoveTriangles(
	const Slicer_2& slicer,
	const std::vector<TRiangle>& remove_triangles) const
{
	std::vector<TRiangle> surface_triangles;
	surface_triangles.reserve(remove_triangles.size());

	const size_t face_cnt = Normals.rows();
	for (size_t i = 0; i < remove_triangles.size(); i++) {
		Eigen::Vector3d v1, v2, v3;
		v1.x() = slicer.positions[remove_triangles[i][0]][0];
		v1.y() = slicer.positions[remove_triangles[i][0]][1];
		v1.z() = slicer.positions[remove_triangles[i][0]][2];
		v2.x() = slicer.positions[remove_triangles[i][1]][0];
		v2.y() = slicer.positions[remove_triangles[i][1]][1];
		v2.z() = slicer.positions[remove_triangles[i][1]][2];
		v3.x() = slicer.positions[remove_triangles[i][2]][0];
		v3.y() = slicer.positions[remove_triangles[i][2]][1];
		v3.z() = slicer.positions[remove_triangles[i][2]][2];
		double ans = (v2.x() - v1.x()) * (v2.y() - v3.y()) - (v2.y() - v1.y()) * (v2.x() - v3.x());
		if (ans > 0) {
			swap(v2, v3);
		}
		double na = (v2.y() - v1.y()) * (v3.z() - v1.z()) - (v2.z() - v1.z()) * (v3.y() - v1.y());
		double nb = (v2.z() - v1.z()) * (v3.x() - v1.x()) - (v2.x() - v1.x()) * (v3.z() - v1.z());
		double nc = (v2.x() - v1.x()) * (v3.y() - v1.y()) - (v2.y() - v1.y()) * (v3.x() - v1.x());

		Eigen::Vector3d normal_vector(na, nb, nc);
		normal_vector.normalize();

		bool jud_surface = false;
		for (size_t j = 0; j < face_cnt; j++) {
			double n_dot = normal_vector.dot(Normals.row(j));
			if (fabs(n_dot) < 0.999) {
				continue;
			}
			Eigen::Vector3d v1_s(v1.x(), v1.y(), v1.z());
			Eigen::Vector3d v2_s(v2.x(), v2.y(), v2.z());
			Eigen::Vector3d v3_s(v3.x(), v3.y(), v3.z());

			bool b1 = vasco::core::isPointInTriangle(v1_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));
			bool b2 = vasco::core::isPointInTriangle(v2_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));
			bool b3 = vasco::core::isPointInTriangle(v3_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));

			if (b1 && b2 && b3) {
				jud_surface = true;
				break;
			}
		}

		if (jud_surface) {
			surface_triangles.push_back(remove_triangles[i]);
		}
	}

	return surface_triangles;
}

void HybridManufacturing::InitializeVoronoi()
{

	std::vector<vasco::VoronoiCell> voronoiCells;
	std::vector<Eigen::Vector3d> bottomVertices;

	std::vector<vasco::VoronoiCell> voronoiCells1;
	std::vector<Eigen::Vector3d> bottomVertices1;
	vasco::BuildVoronoiCells(V, F, 2.0, voronoiCells, bottomVertices,
		open_vis_voronoi, file_name);

	vasco::BuildVoronoiCells(input_mesh, 2.0, voronoiCells1, bottomVertices1,
		open_vis_voronoi, file_name);

	size_t cell_cnt = voronoiCells.size();
	for (size_t i = 0; i < cell_cnt; i++) {
		if (voronoiCells[i].is_available != voronoiCells1[i].is_available) {
			std::cerr << "Error: Voronoi cell availability mismatch at index " << i << std::endl;
		}
		if (voronoiCells[i].is_available) {
			size_t point_cnt = voronoiCells[i].all_points_in_polygon.size();
			if (point_cnt != voronoiCells1[i].all_points_in_polygon.size()) {
				std::cerr << "Error: Voronoi cell point count mismatch at index " << i << std::endl;
			}
			//for (size_t j = 0; j < point_cnt; j++) {
			//	if ((voronoiCells[i].all_points_in_polygon[j] - voronoiCells1[i].all_points_in_polygon[j]).norm() > 1e-6) {
			//		std::cerr << "Error: Voronoi cell point mismatch at index " << i << ", point " << j << std::endl;
			//	}
			//}
			if (voronoiCells[i].site != voronoiCells1[i].site) {
				std::cerr << "Error: Voronoi cell site mismatch at index " << i << std::endl;
			}
			if (voronoiCells[i].adjacent_cells.size() != voronoiCells1[i].adjacent_cells.size()) {
				std::cerr << "Error: Voronoi cell adjacent cell count mismatch at index " << i << std::endl;
			}
		}
	}
	all_voronoi_cells = voronoiCells;
	V_bottom = bottomVertices;
}

void HybridManufacturing::InitializePolyscope()
{
	polyscope::init();
	polyscope::view::setUpDir(polyscope::UpDir::ZUp);
}

int HybridManufacturing::CollisionDetectionForSubtractiveManufacturing(cutter the_nozzle)
{
	clock_t start_time, end_time;
	clock_t start_time_2, end_time_2;
	clock_t start_time_3, end_time_3;
	clock_t start_time_4, end_time_4;
	double sum_time = 0;
	int test_cont_num = 0;

	start_time = clock();
	std::cout << "Doing collision detection for subtractive manufacturing......" << endl;
	string file_name_2 = file_name;
	int cont_num = 0;
	sampling_subtractive.OrientationSamplePoints();	//sampling_subtractive生成球面采样点
	Eigen::Matrix3d rotMatrix;
	std::vector<bool> flag_accessible_points;	//true值代表对应的voronoi单元在某个采样方向下可达
	std::vector<bool> flag_covering_points;	//true值代表对应的voronoi单元在某个采样方向下被喷嘴覆盖
	//double nozzle_par = (the_nozzle.upper_surface_r - the_nozzle.lowwer_surface_r) / the_nozzle.nozzle_H_half;

	flag_accessible_points.resize(all_voronoi_cells.size());
	flag_covering_points.resize(all_voronoi_cells.size());
	flag_voronoi_has_been_printed.resize(all_voronoi_cells.size());
	for (int i = 0; i < all_voronoi_cells.size(); i++) {
		flag_accessible_points[i] = false;
		flag_voronoi_has_been_printed[i] = true;
		flag_covering_points[i] = false;
	}

	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {	//枚举所有采样方向

		///////////////////rotate/////////////////////
		std::vector<Eigen::Vector3d> temp_V;
		temp_V.resize(V.rows());	//temp_V存储旋转后的模型顶点
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i].x() = V.row(i).x();
			temp_V[i].y() = V.row(i).y();
			temp_V[i].z() = V.row(i).z();
		}
		vector<std::vector<Eigen::Vector3d>> temp_new_V;	//temp_new_V存储旋转后的voronoi面边界顶点
		temp_new_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_new_V[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
			for (int k = 0; k < temp_new_V[i].size(); k++) {
				temp_new_V[i][k].x() = all_voronoi_cells[i].all_points_in_polygon[k].x();
				temp_new_V[i][k].y() = all_voronoi_cells[i].all_points_in_polygon[k].y();
				temp_new_V[i][k].z() = all_voronoi_cells[i].all_points_in_polygon[k].z();
			}
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();	//计算从z轴到采样方向的旋转矩阵
		Eigen::Matrix3d rotMatrix_inverse = rotMatrix.inverse();
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix_inverse * temp_V[i];
		for (int i = 0; i < V.rows(); i++)
			for (int j = 0; j < temp_new_V[i].size(); j++)
				temp_new_V[i][j] = rotMatrix_inverse * temp_new_V[i][j];
		//////////////////////////////////////////////


		/////////////////calculate normal//////////////////
		//vector<vector<Vector3>> all_normal_of_triangles_in_cells;
		vector<Eigen::Vector3d> all_normal_of_cells;	//存储每个voronoi单元的法向量
		//all_normal_of_triangles_in_cells.resize(all_voronoi_cells.size());
		all_normal_of_cells.resize(all_voronoi_cells.size());

		vector<vector<Eigen::Vector3d>> temp_vis;
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			const auto& cell = all_voronoi_cells[i];
			if (cell.is_available == true) {
				const auto polygon_size = cell.all_points_in_polygon.size();
				const int site_index = cell.site;

				Eigen::Vector3d vn_sum(0.0, 0.0, 0.0);
				for (int j = 0; j < polygon_size; j++) {
					const Eigen::Vector3d v1 = ToVector3(temp_V[site_index]);
					const Eigen::Vector3d v2 = ToVector3(temp_new_V[i][j]);
					const Eigen::Vector3d v3 = ToVector3(temp_new_V[i][(j + 1) % polygon_size]);
					vn_sum += (v2 - v1).cross(v3 - v1);
				}

				all_normal_of_cells[i] = vn_sum;
				if (polygon_size > 0) {
					all_normal_of_cells[i] /= polygon_size;
				}
				all_normal_of_cells[i].normalize();	//voronoi单元法向量归一化

				all_voronoi_cells[i].all_normal_in_all_ori.push_back(all_normal_of_cells[i]);	//将该采样方向下的法向量存入all_voronoi_cells的all_normal_in_all_ori中

				temp_vis.push_back({});
				temp_vis.back().reserve(2);
				Eigen::Vector3d v_site(
					temp_V[site_index].x(),
					temp_V[site_index].y(),
					temp_V[site_index].z());
				Eigen::Vector3d v_normal = v_site + all_normal_of_cells[i] * 3.0;
				temp_vis.back().push_back(v_site);
				temp_vis.back().push_back(v_normal);
			}
		}
		if (ori == 0) {
			Visual vis_normal;//存储法向量的可视化
			//vis_normal.generateModelForRendering_9(temp_vis, file_name);

		}
		//////////////////////////////////////////////////


		///////////////////collision detection////////////////////////
		int cont_accessible_points = 0, cont_unaccessible_points = 0;

		vector<double> max_z_of_cells(all_voronoi_cells.size());	//max_z_of_cells存储每个voronoi单元边界顶点(all_points_in_polygon)的最大z值
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			max_z_of_cells[i] = MIN_D;
			for (int j = 0; j < all_voronoi_cells[i].all_points_in_polygon.size(); j++)
				max_z_of_cells[i] = max(max_z_of_cells[i], temp_new_V[i][j](2, 0));
		}

		PrepareToolForCollision(the_nozzle);

		for (int i = 0; i < all_voronoi_cells.size(); i++) {	//枚举所有voronoi单元，判断该单元在当前采样方向ori下是否可达
			if (flag_accessible_points[i] == true) {	//如果该单元在之前某个采样方向下已经可达，则跳过
				cont_accessible_points++;
				continue;
			}
			cont_unaccessible_points++;

			if (!all_voronoi_cells[i].is_available) {
				continue;
			}

			// 计算刀尖球中心点坐标
			Eigen::Vector3d center_point = ComputeToolCenter(
				ToVector3(temp_V[all_voronoi_cells[i].site]),
				all_normal_of_cells[i],
				the_nozzle.cylinder_r);

			bool jud_collision = false;
			for (int ii = 0; ii < all_voronoi_cells.size(); ii++) {
				if (i == ii || !all_voronoi_cells[ii].is_available) {
					continue;
				}

				if (CheckToolCollisionWithCell(center_point, temp_new_V[ii], max_z_of_cells[ii], the_nozzle)) {
					jud_collision = true;
					break;
				}
			}

			if (!jud_collision) {
				flag_accessible_points[i] = true;
			}
		}
		num_inaccessible_points = cont_unaccessible_points;
		std::cout << "id of orientation:" << ori << endl;
		std::cout << "number of accessible points:" << cont_accessible_points << endl;
		std::cout << "number of unaccessible points:" << cont_unaccessible_points << endl << endl;
		//break;
	}
	end_time = clock();


	//*******************************************//
	/////////////////find area S//////////////////
	//*******************************************//
	start_time_2 = clock();

	int cont_number_2 = 0;
	std::vector<Eigen::MatrixXd> vis_red_points;	//vis_red_points存储不可达voronoi单元的边界顶点坐标，用于可视化
	std::vector<Eigen::MatrixXd> vis_green_points;	//vis_green_points存储可达voronoi单元的边界顶点坐标，用于可视化
	std::vector<Eigen::Vector3d> temp_V_vis;	//temp_V_vis存储未旋转的模型顶点坐标，用于可视化
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {	//ori枚举所有采样方向
		std::vector<Eigen::Vector3d> temp_V;	//temp_V存储旋转后的模型顶点坐标
		temp_V.resize(V.rows());
		temp_V_vis.resize(V.rows());	//temp_V_vis存储未旋转的模型顶点坐标
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i](0, 0) = V.row(i).x();
			temp_V[i](1, 0) = V.row(i).y();
			temp_V[i](2, 0) = V.row(i).z();
		}
		vector<std::vector<Eigen::Vector3d>> temp_new_V;	//temp_new_V存储旋转后的voronoi面边界顶点坐标
		temp_new_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_new_V[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
			for (int k = 0; k < temp_new_V[i].size(); k++) {
				temp_new_V[i][k](0, 0) = all_voronoi_cells[i].all_points_in_polygon[k].x();
				temp_new_V[i][k](1, 0) = all_voronoi_cells[i].all_points_in_polygon[k].y();
				temp_new_V[i][k](2, 0) = all_voronoi_cells[i].all_points_in_polygon[k].z();
			}
		}
		temp_V_vis = temp_V;

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix.inverse() * temp_V[i];	//将模型顶点旋转到以采样方向为z轴的坐标系下
		for (int i = 0; i < V.rows(); i++)
			for (int j = 0; j < temp_new_V[i].size(); j++)
				temp_new_V[i][j] = rotMatrix.inverse() * temp_new_V[i][j];	//将voronoi面边界顶点旋转到以采样方向为z轴的坐标系下

		//cout << "b" << endl;
		/////////////////calculate normal//////////////////
		//vector<vector<Vector3>> all_normal_of_triangles_in_cells;
		vector<Eigen::Vector3d> all_normal_of_cells;	//仅计算不可达voronoi单元的法向量
		//all_normal_of_triangles_in_cells.resize(all_voronoi_cells.size());
		all_normal_of_cells.resize(all_voronoi_cells.size());
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			if (all_voronoi_cells[i].is_available == true && flag_accessible_points[i] == false) {	//仅计算不可达voronoi单元的法向量
				//all_normal_of_triangles_in_cells[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
				all_normal_of_cells[i].x() = all_normal_of_cells[i].y() = all_normal_of_cells[i].z() = 0;
				for (int j = 0; j < 1; j++) {
					Eigen::Vector3d v1 = ToVector3(temp_V[all_voronoi_cells[i].site]);
					Eigen::Vector3d v2 = ToVector3(temp_new_V[i][j]);
					Eigen::Vector3d v3 = ToVector3(temp_new_V[i][(j + 1) % all_voronoi_cells[i].all_points_in_polygon.size()]);
					Eigen::Vector3d vn = (v2 - v1).cross(v3 - v1);
					//all_normal_of_triangles_in_cells[i][j] = vn;
					all_normal_of_cells[i].x() += vn.x();
					all_normal_of_cells[i].y() += vn.y();
					all_normal_of_cells[i].z() += vn.z();
				}
				all_normal_of_cells[i] /= all_voronoi_cells[i].all_points_in_polygon.size();
				all_normal_of_cells[i].normalize();
			}
		}
		//////////////////////////////////////////////////
		//cout << "c" << endl;
		//ofstream vis_unaccessible_points(file_name + "_unaccessible_points.obj");

		file_name_2 = file_name + to_string(cont_num);
		vector<double> max_z_of_cells(all_voronoi_cells.size());	//max_z_of_cells存储每个voronoi单元边界顶点(all_points_in_polygon)的最大z值
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			max_z_of_cells[i] = MIN_D;
			for (int j = 0; j < all_voronoi_cells[i].all_points_in_polygon.size(); j++)
				max_z_of_cells[i] = max(max_z_of_cells[i], temp_new_V[i][j](2, 0));
		}

		int cont_unaccessible_points = 0;

		PrepareToolForCollision(the_nozzle);

		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			if (flag_accessible_points[i] == true) {	//跳过可达voronoi单元
				continue;
			}

			if (!all_voronoi_cells[i].is_available) {
				continue;
			}

			// 计算不可达voronoi单元i的刀尖球中心点坐标
			Eigen::Vector3d center_point = ComputeToolCenter(
				ToVector3(temp_V[all_voronoi_cells[i].site]),
				all_normal_of_cells[i],
				the_nozzle.cylinder_r);

			int index_insert = 0;
			if (ori == 0) {	//在首个方向 ori==0 时，初始化该不可达单元的 area_S 容器与索引映射
				vector<area_S> temp_vec_area_s;
				//temp_vec_area_s.resize(1);
				all_the_area_S.push_back(temp_vec_area_s);
				map_S_and_vertex.insert({ cont_unaccessible_points,i });
				map_S_and_vertex_inv.insert({ i, cont_unaccessible_points });
				cont_unaccessible_points++;	//记录不可达单元数量,第(cont_unaccessible_points-1)个不可达voronoi单元的索引为i
			}
			else {
				index_insert = map_S_and_vertex_inv[i];
			}

			//vis_unaccessible_points << "v " << temp_V[i](0,0) << " " << temp_V[i](1, 0) << " " << temp_V[i](2, 0) << endl;
			start_time_4 = clock();
			if (ori == 0)
				vis_red_points.push_back(temp_V_vis[i]);	//将不可达voronoi单元i的网格顶点存入vis_red_points，用于可视化
			end_time_4 = clock();
			sum_time += (double(end_time_4) - double(start_time_4)) / CLOCKS_PER_SEC;

			for (int ii = 0; ii < all_voronoi_cells.size(); ii++) {	//枚举所有voronoi单元，判断刀尖放在当前单元i、刀具方向是ori方向时，是否与其他单元发生碰撞
				if (i == ii || !all_voronoi_cells[ii].is_available) {
					continue;
				}

				if (CheckToolCollisionWithCell(center_point, temp_new_V[ii], max_z_of_cells[ii], the_nozzle)) {
					area_S temp_area_S(ii, ori);
					if (ori == 0) {
						all_the_area_S.back().push_back(temp_area_S);
					}
					else {
						all_the_area_S[index_insert].push_back(temp_area_S);
						cont_number_2++;
					}
					test_cont_num++;
				}
			}
		}

		//cout << "d" << endl;
		Visual Vis;
		//Vis.generateModelForRendering(temp_layers, file_name_2);
		//Vis.generateModelForRendering_3(vectorAfter, file_name_2, vis_points);
		//break;

		cont_num++;
	}
	//*******************************************//
	//cout << endl << all_the_area_S.size() << endl;
	//cout << endl << cont_number_2 << endl;
	end_time_2 = clock();

	////////////get covering points(green points)////////////
	//only consider covering for area S,save the index of all_the_area_S
	start_time_3 = clock();
	cout << test_cont_num << "&&&&" << endl;
	int cont_covering_points = 0;
	for (int i = 0; i < all_the_area_S.size(); i++) {	//枚举第i个不可达voronoi单元的area_S
		for (int j = 0; j < all_the_area_S[i].size(); j++) {	//枚举里面的所有area_S元素
			if (flag_covering_points[all_the_area_S[i][j].pointId] == false) {	//如果该voronoi单元还没有在flag_covering_points中被标记
				map_covering_points_and_vertex.insert({ cont_covering_points,all_the_area_S[i][j].pointId });	//建立被碰撞点与voronoi单元索引的映射
				map_covering_points_and_vertex_inv.insert({ all_the_area_S[i][j].pointId, cont_covering_points });
				vector<area_S> temp_vec_area_s;
				all_the_covering_points.push_back(temp_vec_area_s);	//在all_the_covering_points中添加一个新的容器

				area_S temp_covering_point(i, all_the_area_S[i][j].oriId);	//
				all_the_covering_points[all_the_covering_points.size() - 1].push_back(temp_covering_point);
				flag_covering_points[all_the_area_S[i][j].pointId] = true;
				cont_covering_points++;
			}
			else {
				int index_insert = map_covering_points_and_vertex_inv[all_the_area_S[i][j].pointId];
				area_S temp_covering_point(i, all_the_area_S[i][j].oriId);
				all_the_covering_points[index_insert].push_back(temp_covering_point);
			}
		}
	}
	int max_size = -MAX_I;
	for (int i = 0; i < all_the_covering_points.size(); i++) {
		max_size = max(max_size, int(all_the_covering_points[i].size()));
		vis_green_points.push_back(temp_V_vis[map_covering_points_and_vertex[i]]);	//将覆盖点对应的voronoi单元顶点存入vis_green_points，用于可视化
	}
	vector<double> color_map;
	color_map.resize(all_the_covering_points.size());
	for (int i = 0; i < all_the_covering_points.size(); i++)
		color_map[i] = double(all_the_covering_points[i].size()) / double(max_size);	//归一化颜色值，便于可视化，被更多个area S覆盖的点颜色更亮
	end_time_3 = clock();
	/*vector<vector<area_S>> temp_all_the_covering_points = all_the_covering_points;
	map<int, int> temp_map_covering_points_and_vertex = map_covering_points_and_vertex;
	for(int i = 0;i< temp_all_the_covering_points.size();i++)
		for (int j = i + 1; j < temp_all_the_covering_points.size(); j++)
			if (temp_all_the_covering_points[i].size() < temp_all_the_covering_points[j].size()) {
				swap(temp_all_the_covering_points[i], temp_all_the_covering_points[j]);
				swap(temp_map_covering_points_and_vertex[i], temp_map_covering_points_and_vertex[j]);
			}

	for(int i =0;i < temp_all_the_covering_points.size() / 10;i++)
		vis_green_points.push_back(temp_V_vis[temp_map_covering_points_and_vertex[i]]);*/
		/////////////////////////////////////////////////////////
	if (open_vis_red_points == true)
		createRedBalls(file_name, vis_red_points);
	if (open_vis_green_points == true)
		createGreenBalls(file_name, vis_green_points, color_map);

	std::cout << "Collision Detection For Subtractive Manufacturing have done" << endl;
	/*std::cout << "&&&time&&& Collision detection: " << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
	std::cout << "&&&time&&& Subtractive dependency graph: " << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << std::endl;
	std::cout << "&&&time&&& Green points generation: " << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << std::endl;*/
	time_build_subtractive_graph = double(end_time - start_time) / CLOCKS_PER_SEC + double(end_time_2 - start_time_2) / CLOCKS_PER_SEC + double(end_time_3 - start_time_3) / CLOCKS_PER_SEC;
	std::cout << "&&&time&&& Build subtractive graph: " << double(end_time - start_time) / CLOCKS_PER_SEC + double(end_time_2 - start_time_2) / CLOCKS_PER_SEC + double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << std::endl;

	if (all_the_area_S.size() == 0) {
		cout << file_name << " 无不可达点！！！" << endl;

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXi N;

		igl::readOBJ(file_name + ".obj", V, F);
		igl::writeOBJ(file_name + "-1_0_current.obj", V, F);
		return -1;
	}
	return 0;
}

void HybridManufacturing::GetALLFragileVertex(SAMPLE_ON_BALL sampling)
{
	for (int ori = 0; ori < sampling.sample_points.size(); ori++) {
		//rotating the blocks and then slicing//
		vector<Eigen::Vector3d> temp_V;
		temp_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i](0, 0) = V.row(i).x();
			temp_V[i](1, 0) = V.row(i).y();
			temp_V[i](2, 0) = V.row(i).z();
		}
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling.sample_points[ori]);
		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix.inverse() * temp_V[i];

		vector<Eigen::Vector3d> normal_V;
		vector<bool> temp_fragile_V;
		normal_V.resize(V.rows());
		temp_fragile_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			normal_V[i].x() = normal_V[i].y() = normal_V[i].z() = 0;
		}
		for (int i = 0; i < F.rows(); i++) {
			VEctor v1, v2, v3;
			v1[0] = temp_V[F(i, 0)](0, 0);
			v1[1] = temp_V[F(i, 0)](1, 0);
			v1[2] = temp_V[F(i, 0)](2, 0);
			v2[0] = temp_V[F(i, 1)](0, 0);
			v2[1] = temp_V[F(i, 1)](1, 0);
			v2[2] = temp_V[F(i, 1)](2, 0);
			v3[0] = temp_V[F(i, 2)](0, 0);
			v3[1] = temp_V[F(i, 2)](1, 0);
			v3[2] = temp_V[F(i, 2)](2, 0);
			Eigen::Vector3d temp_normal;
			temp_normal = faceNormal(v1, v2, v3);
			temp_normal.normalize();
			normal_V[F(i, 0)] += temp_normal;
			normal_V[F(i, 1)] += temp_normal;
			normal_V[F(i, 2)] += temp_normal;
		}
		for (int i = 0; i < V.rows(); i++) {
			normal_V[i].normalize();
			if (normal_V[i].z() > 0.95)
				temp_fragile_V[i] = true;
			else
				temp_fragile_V[i] = false;
		}
		is_fragile_V.push_back(temp_fragile_V);
	}
}

void HybridManufacturing::detect_collision_with_printing_platform(
	int& index,
	vector<int>& candidate_nodes,
	OrientationScores& all_calculated_value,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d ori_now,
	nozzle the_nozzle)
{

	Eigen::Vector3d center_point(0.0, 0.0, 0.0);
	for (int i = 0; i < V_bottom.size(); i++) {
		center_point += V_bottom[i];
	}
	center_point /= V_bottom.size();
	double circle_r = 160;

	if (ori_now.x() == 0 && ori_now.y() == 0 && ori_now.z() == 1)
		return;
	//if (all_calculated_value[index].number_of_remaining_face < terminate_threshold_of_number_of_faces) 
		//return;
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(ori_now, vectorBefore).toRotationMatrix();
	vector<vector<Eigen::Vector3d>> rotate_all_cut_layers;
	rotate_all_cut_layers.resize(all_cut_layers.size());
	for (int i = 0; i < all_cut_layers.size(); i++) {
		rotate_all_cut_layers[i].resize(all_cut_layers[i].size());
		for (int j = 0; j < all_cut_layers[i].size(); j++) {
			rotate_all_cut_layers[i][j](0, 0) = all_cut_layers[i][j].x;
			rotate_all_cut_layers[i][j](1, 0) = all_cut_layers[i][j].y;
			rotate_all_cut_layers[i][j](2, 0) = all_cut_layers[i][j].z;
			rotate_all_cut_layers[i][j] = rotMatrix.inverse() * rotate_all_cut_layers[i][j];
		}
	}

	for (int i = 0; i < rotate_all_cut_layers.size(); i++) {
		for (int j = 0; j < rotate_all_cut_layers[i].size(); j++) {
			if ((abs(rotate_all_cut_layers[i][j](2, 0) - center_point.z()) < the_nozzle.lowwer_surface_r) && (pow(rotate_all_cut_layers[i][j](0, 0) - center_point.x(), 2) + pow(rotate_all_cut_layers[i][j](1, 0) - center_point.y(), 2) - pow(circle_r, 2) < 0)) {
				candidate_nodes.erase(candidate_nodes.begin() + index);
				all_calculated_value.erase(all_calculated_value.begin() + index);
				index--;
				return;
			}
		}
	}
}

all_value HybridManufacturing::GainMesh(
	Slicer_2& slicer,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d vector_after,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	vector<int> all_cut_layers_dependency_layer,
	bool flag_is_continue_block,
	int id_continue)
{
	clock_t start_time, end_time;

	double sum_area = 0;
	double sum_value_of_manufacturing_require = 0;
	double sum_not_self_support_area = 0;
	//double area_threshold = 5.0;
	double ratio_of_area_threshold = 0.05;  //0.015 //0.15  //0.08   //0.03  //0.05
	double max_self_support_slope_angle;
	//cout << "#####" << height_of_beam_search << " " << id_continue << endl;
	//if(height_of_beam_search != 7 || id_continue != 1)s
	//	max_self_support_slope_angle = PI / 10;  //PI/3.6  
	//else
	max_self_support_slope_angle = PI / 3.6;  //PI/3.6  
	double sum_area_threshold = 1000;  //200 //1000
	if (vector_after.x() == 0 && vector_after.y() == 0) {
		sum_area_threshold = 100000;
		ratio_of_area_threshold = 0.1;
	}

	/*if(height_of_beam_search == 4)
		ratio_of_area_threshold = 0.2;*/

	std::vector<bool> jud_triangle_have_been_added;
	vector<vector<int>> cutting_plane_points;
	cutting_plane_points.resize(all_cut_layers.size());
	//////sort cut layers//////
	for (int i = 0; i < all_cut_layers.size(); i++)
		for (int j = i + 1; j < all_cut_layers.size(); j++) {

			if (all_cut_layers[i][0].z > all_cut_layers[j][0].z) {
				swap(all_cut_layers[i], all_cut_layers[j]);
				swap(all_cut_layers_dependency_layer[i], all_cut_layers_dependency_layer[j]);
			}
		}
	///////////////////////////


	//need rotate first//
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	for (int i = 0; i < slicer.positions.size(); i++) {
		Eigen::Vector3d temp_V;
		temp_V(0, 0) = slicer.positions[i][0];
		temp_V(1, 0) = slicer.positions[i][1];
		temp_V(2, 0) = slicer.positions[i][2];
		temp_V = rotMatrix.inverse() * temp_V;
		slicer.positions[i][0] = temp_V(0, 0);
		slicer.positions[i][1] = temp_V(1, 0);
		slicer.positions[i][2] = temp_V(2, 0);
	}
	/////////////////////


	//--------------------cut-----------------------//
	slicer.normal[0] = 0;
	slicer.normal[1] = 0;
	slicer.normal[2] = 1;
	//cout << endl << "pp" << endl;
	for (int i = 0; i < all_cut_layers.size(); i++) {
		slicer.origin[0] = 0;
		slicer.origin[1] = 0;
		slicer.origin[2] = all_cut_layers[i][0].z;
		slicer.cut();
	}
	//cout << endl << "qq" << endl;
	std::vector<TRiangle> ori_triangle = slicer.triangles;
	jud_triangle_have_been_added.resize(slicer.triangles.size());
	for (int i = 0; i < slicer.triangles.size(); i++)
		jud_triangle_have_been_added[i] = false;

	Slicer_2 all_slicer;

	all_slicer.positions = slicer.positions;
	all_slicer.triangles = ori_triangle;


	//--------------------save candidate_triangles-----------------------//
	int current_index = 0;
	std::vector<VEctor> min_z_point;
	std::vector<double> min_z_triangle;
	vector<int> index_of_min_point_in_triangle;
	std::vector<TRiangle> candidate_triangles;
	std::vector<int> id_triangles;
	std::vector<TRiangle> remove_triangles;
	candidate_triangles.clear();
	id_triangles.clear();
	double min_z_all_cut_layers = 9999999;
	for (int t = 0; t < all_cut_layers.size(); t++)
		min_z_all_cut_layers = min(min_z_all_cut_layers, all_cut_layers[t][0].z);
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;
		for (int k = 0; k < 3; k++) {
			temp_min_z_triangle = min(all_slicer.positions[all_slicer.triangles[i][k]][2], temp_min_z_triangle);
		}
		if (temp_min_z_triangle + 0.001 >= min_z_all_cut_layers) {
			candidate_triangles.push_back(all_slicer.triangles[i]);
			id_triangles.push_back(i);
		}
	}
	for (int i = 0; i < candidate_triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;
		VEctor temp_min_z_point = all_slicer.positions[candidate_triangles[i][0]];
		int temp_index_of_min_point_in_triangle;
		for (int k = 0; k < 3; k++) {
			if (all_slicer.positions[candidate_triangles[i][k]][2] < temp_min_z_triangle) {
				temp_min_z_triangle = all_slicer.positions[candidate_triangles[i][k]][2];
				temp_min_z_point = all_slicer.positions[candidate_triangles[i][k]];
				temp_index_of_min_point_in_triangle = candidate_triangles[i][k];
			}
		}
		min_z_triangle.push_back(temp_min_z_triangle);
		min_z_point.push_back(temp_min_z_point);
		index_of_min_point_in_triangle.push_back(temp_index_of_min_point_in_triangle);
	}

	quick_sort(candidate_triangles, 0, candidate_triangles.size() - 1, id_triangles, min_z_triangle, min_z_point);


	////--------------------save OPP_triangles one by one-----------------------//
	vector<vector<TRiangle>> all_furcation_of_blocks;
	all_furcation_of_blocks.resize(all_cut_layers.size());
	vector<int> id_of_furcation_of_blocks;
	vector<int> save_current_index;
	for (int t = 0; t < all_cut_layers.size(); t++) {
		current_index = 0;   //该方式也许比较慢
		double boundary_bottom = 999999, boundary_left = 999999, boundary_top = -999999, boundary_right = -999999;
		for (int i = 0; i < all_cut_layers[t].size(); i++) {
			boundary_top = std::max(boundary_top, all_cut_layers[t][i].y);
			boundary_bottom = std::min(boundary_bottom, all_cut_layers[t][i].y);
			boundary_right = std::max(boundary_right, all_cut_layers[t][i].x);
			boundary_left = std::min(boundary_left, all_cut_layers[t][i].x);
		}
		Eigen::Vector2d current_triangle_point;
		Eigen::Vector2d current_layer_point;
		for (; current_index < candidate_triangles.size(); current_index++) {
			if (abs(min_z_triangle[current_index] - all_cut_layers[t][0].z) > 0.001) {
				if (min_z_triangle[current_index] > all_cut_layers[t][0].z)
					break;
			}
			else {
				bool jud_is_boundary_point = false;
				//for (int i = 0; i < all_cut_layers[t].size(); i++) {
				int cont_inside_boundary = 0;
				for (int k = 0; k < 3; k++) {
					if (all_slicer.positions[candidate_triangles[current_index][k]][0] + 0.01 >= boundary_left && all_slicer.positions[candidate_triangles[current_index][k]][0] - 0.01 <= boundary_right
						&& all_slicer.positions[candidate_triangles[current_index][k]][1] + 0.01 >= boundary_bottom && all_slicer.positions[candidate_triangles[current_index][k]][1] - 0.01 <= boundary_top) {
						cont_inside_boundary++;
					}
				}
				if (cont_inside_boundary >= 2)
					jud_is_boundary_point = true;
				/*if (min_z_point[current_index][0] + 0.01 >= boundary_left && min_z_point[current_index][0] - 0.01 <= boundary_right &&
					min_z_point[current_index][1] + 0.01 >= boundary_bottom && min_z_point[current_index][1] - 0.01 <= boundary_top) {
					jud_is_boundary_point = true;
					break;
				}*/
				//}
				if (jud_is_boundary_point == true) {
					for (int j = 0; j < all_cut_layers[t].size(); j++) {
						current_layer_point.x() = all_cut_layers[t][j].x;
						current_layer_point.y() = all_cut_layers[t][j].y;
						current_triangle_point.x() = min_z_point[current_index][0];
						current_triangle_point.y() = min_z_point[current_index][1];
						if ((current_layer_point - current_triangle_point).norm() < 0.2) {  //4.0
							//if (height_of_beam_search != 2)
							remove_triangles.push_back(candidate_triangles[current_index]);
							id_of_furcation_of_blocks.push_back(t);
							all_furcation_of_blocks[t].push_back(candidate_triangles[current_index]);

							Eigen::Vector3d face_normal = faceNormal(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							face_normal.normalize();
							Eigen::Vector3d base_normal(0, 0, 1);
							int jud_self_support = (face_normal.dot(base_normal) + sin(max_self_support_slope_angle) >= 0);
							double current_area = triangleArea(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							sum_area += current_area;
							//sum_value_of_manufacturing_require += current_area * calculate_manufacturing_require(all_slicer, remove_triangles);
							if (jud_self_support == false) {
								sum_not_self_support_area += current_area;
								if (sum_not_self_support_area / sum_area > ratio_of_area_threshold && sum_area >= sum_area_threshold) {
									all_value temp_all_value;
									temp_all_value.value_of_self_support = 0;
									return temp_all_value;
								}
							}

							save_current_index.push_back(current_index);
							jud_triangle_have_been_added[id_triangles[current_index]] = true;
							break;
						}
					}
				}
			}

		}
	}
	current_index = 0;  //current_index = 0;
	//std::cout << "**" << current_index << " " << candidate_triangles.size() << " " << remove_triangles.size() << endl;

	start_time = clock();
	while (1) {
		bool flag_break = false;
		for (int i = 0; i < remove_triangles.size(); i++)
		{
			if (jud_triangle_have_been_added[id_triangles[current_index]] == false && min_z_triangle[current_index] >= min_z_triangle[save_current_index[i]]) {
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						if (candidate_triangles[current_index][j] == remove_triangles[i][k]) {
							//if(height_of_beam_search != 2)
							remove_triangles.push_back(candidate_triangles[current_index]);
							id_of_furcation_of_blocks.push_back(id_of_furcation_of_blocks[i]);
							all_furcation_of_blocks[id_of_furcation_of_blocks[i]].push_back(candidate_triangles[current_index]);
							//////all_furcation_of_blocks[t].push_back(candidate_triangles[current_index]);

							Eigen::Vector3d face_normal = faceNormal(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							face_normal.normalize();
							Eigen::Vector3d base_normal(0, 0, 1);
							int jud_self_support = (face_normal.dot(base_normal) + sin(max_self_support_slope_angle) >= 0);
							double current_area = triangleArea(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							sum_area += current_area;

							if (jud_self_support == false) {
								sum_not_self_support_area += current_area;
								if (sum_not_self_support_area / sum_area > ratio_of_area_threshold && sum_area >= sum_area_threshold) {
									all_value temp_all_value;
									temp_all_value.value_of_self_support = 0;
									return temp_all_value;
								}
							}

							save_current_index.push_back(current_index);
							jud_triangle_have_been_added[id_triangles[current_index]] = true;
							flag_break = true;
							break;
						}
					}
					if (flag_break == true)
						break;
				}
			}
			if (flag_break == true)
				break;
		}
		current_index++;
		if (current_index >= candidate_triangles.size())
			break;
	}end_time = clock();
	//add remaining face, set a distance threshold Dis, only a face less Dis from other cut layers and no other dependency layer exist, the face is consider to remaining face
	/*double Dis = dh * 4;
	for (int i = 0; i < jud_triangle_have_been_added.size(); i++) {
		if (jud_triangle_have_been_added[i] == false) {
			int id_layer;
			double min_dis = 99999999;
			for (int t = 0; t < all_cut_layers.size(); t++) {
				if (all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] <= 2 * dh && all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] >= -0.001) {
					for (int j = 0; j < all_cut_layers[t].size(); j++) {
						cv::Point3d current_triangle_point(slicer.positions[slicer.triangles[i][0]][0], slicer.positions[slicer.triangles[i][0]][1], slicer.positions[slicer.triangles[i][0]][2]);
						cv::Point3d current_layer_point(all_cut_layers[t][j].x, all_cut_layers[t][j].y, all_cut_layers[t][0].z);

						double distance = distance_3d(current_triangle_point, current_layer_point);
						if (distance < min_dis) {
							min_dis = distance;
							id_layer = t;
						}
					}
				}
			}
			if (min_dis < Dis && all_cut_layers_dependency_layer[id_layer] == 0)
				remove_triangles.push_back(slicer.triangles[i]);
		}
	}*/

	all_value all_calculated_value;
	all_calculated_value.number_of_remaining_face = all_slicer.triangles.size() - remove_triangles.size();
	all_slicer.triangles = remove_triangles;


	all_calculated_value.value_of_projected = calculate_projected_area(all_slicer, all_furcation_of_blocks, all_cut_layers);

	double min_z = 9999999, max_z = -9999999;
	double value_of_manufacturing_require;
	for (int i = 0; i < all_slicer.positions.size(); i++) {
		min_z = min(min_z, all_slicer.positions[i][2]);
		max_z = max(max_z, all_slicer.positions[i][2]);
	}
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		double current_area = triangleArea(all_slicer.positions[all_slicer.triangles[i][0]], all_slicer.positions[all_slicer.triangles[i][1]], all_slicer.positions[all_slicer.triangles[i][2]]);
		value_of_manufacturing_require = vasco::calculateManufacturingRequire(all_slicer, i, min_z, max_z);
		sum_value_of_manufacturing_require += value_of_manufacturing_require * current_area;
	}

	all_calculated_value.value_of_self_support = 1 - sum_not_self_support_area / sum_area;
	all_calculated_value.large_base = sum_value_of_manufacturing_require / sum_area;


	//将slicer旋转回原始位置
	RotateSlicerPositions(all_slicer, vector_after, vectorBefore);

	//auto filtered_patch = FilterSurfaceRemoveTriangles(all_slicer, remove_triangles);
	//if (!all_slicer.positions.empty() && !filtered_patch.empty()) {
	//	const std::string mesh_name = "all_slicer_" + std::to_string(height_of_beam_search) + "_" + std::to_string(cont_number_of_queue);
	//	polyscope::registerSurfaceMesh(mesh_name, all_slicer.positions, filtered_patch);
	//}

	//if (height_of_beam_search > 1) {
	//	subtractive_remove_output(filtered_patch, all_slicer, height_of_beam_search);
	//	cutter cutter1;
	//	int sub_patches = subtractive_accessibility_decomposition_local(height_of_beam_search, cutter1);
	//	std::cout << "[GainMesh] sub_patches: " << sub_patches <<
	//		" cont_number_of_queue " << cont_number_of_queue <<
	//		" index_of_pre_node " << index_of_pre_node <<
	//		" height_of_beam_search " << height_of_beam_search << endl;
	//	all_calculated_value.value_of_sub_patches = 10.0 / sub_patches;
	//	//polyscope::show();

	//	//cout << "()()()(" << double(end_time - start_time) / CLOCKS_PER_SEC << endl;
	//}
	return all_calculated_value;
}

void HybridManufacturing::SortCutLayersByHeight(
	vector<vector<cv::Point3d>>& all_cut_layers,
	vector<int>& all_cut_layers_dependency_layer,
	vector<int>& flag_cut_layers_is_hole,
	map<int, int>& follow_index) const
{
	follow_index.clear();
	for (int i = 0; i < all_cut_layers.size(); i++)
		follow_index.insert({ i,i });	//初始化follow_index
	for (int i = 0; i < all_cut_layers.size(); i++) {
		for (int j = i + 1; j < all_cut_layers.size(); j++) {	//按照z值从小到大排序cut layers
			if (all_cut_layers[i][0].z > all_cut_layers[j][0].z) {
				swap(all_cut_layers[i], all_cut_layers[j]);	//交换cut layers
				swap(follow_index[i], follow_index[j]);	//交换对应关系
				swap(flag_cut_layers_is_hole[i], flag_cut_layers_is_hole[j]);	//交换是否为孔洞的标记
				swap(all_cut_layers_dependency_layer[i], all_cut_layers_dependency_layer[j]);	//交换对应的依赖层数量
			}
		}
	}
}

Slicer_2 HybridManufacturing::LoadSlicerForCutMesh(
	bool flag_is_continue_block,
	int height_of_beam_search,
	int index_of_pre_node,
	int id_continue) const
{
	Slicer_2 slicer;
	if (!flag_is_continue_block) {
		slicer.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + ".obj");
		cout << "&" << endl;
	}
	else {
		slicer.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + "_" + to_string(id_continue - 1) + "_subblock.obj");
		cout << "*" << endl;
	}
	return slicer;
}

void HybridManufacturing::RotateSlicerPositions(
	Slicer_2& slicer,
	const Eigen::Vector3d& vector_before,
	const Eigen::Vector3d& vector_after) const
{
	const Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_before, vector_after).toRotationMatrix();
	for (int i = 0; i < slicer.positions.size(); i++) {
		Eigen::MatrixXd temp_V(3, 1);
		temp_V(0, 0) = slicer.positions[i][0];
		temp_V(1, 0) = slicer.positions[i][1];
		temp_V(2, 0) = slicer.positions[i][2];
		temp_V = rotMatrix.inverse() * temp_V;
		slicer.positions[i][0] = temp_V(0, 0);
		slicer.positions[i][1] = temp_V(1, 0);
		slicer.positions[i][2] = temp_V(2, 0);
	}
}

void HybridManufacturing::RotateLayersForVisualization(
	vector<vector<cv::Point3d>>& all_layers,
	vector<vector<cv::Point3d>>& all_layers_contain,
	const Eigen::Vector3d& vector_after,
	const Eigen::Vector3d& vector_before) const
{
	const Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_after, vector_before).toRotationMatrix();
	for (int i = 0; i < all_layers.size(); i++) {
		for (int j = 0; j < all_layers[i].size(); j++) {
			Eigen::MatrixXd temp_V(3, 1);
			temp_V(0, 0) = all_layers[i][j].x;
			temp_V(1, 0) = all_layers[i][j].y;
			temp_V(2, 0) = all_layers[i][j].z;
			temp_V = rotMatrix.inverse() * temp_V;
			all_layers[i][j].x = temp_V(0, 0);
			all_layers[i][j].y = temp_V(1, 0);
			all_layers[i][j].z = temp_V(2, 0);
		}

		for (int j = 0; j < all_layers_contain[i].size(); j++) {
			Eigen::MatrixXd temp_V(3, 1);
			temp_V(0, 0) = all_layers_contain[i][j].x;
			temp_V(1, 0) = all_layers_contain[i][j].y;
			temp_V(2, 0) = all_layers_contain[i][j].z;
			temp_V = rotMatrix.inverse() * temp_V;
			all_layers_contain[i][j].x = temp_V(0, 0);
			all_layers_contain[i][j].y = temp_V(1, 0);
			all_layers_contain[i][j].z = temp_V(2, 0);
		}
	}
}

void HybridManufacturing::VisualizeCutLayers(
	const vector<vector<cv::Point3d>>& all_layers,
	const vector<vector<cv::Point3d>>& all_layers_contain,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	bool judge_continue_additive,
	int id_continue,
	const Eigen::Vector3d& vector_after) const
{
	Visual Vis;
	cout << id_continue << endl;
	Vis.generateModelForRendering_5(all_layers, all_layers_contain, height_of_beam_search, cont_number_of_queue, file_name, index_of_pre_node, judge_continue_additive, id_continue);
	string vis_file(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "(" + to_string(index_of_pre_node) + ")_Layer");
	string vis_file_contain(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "(" + to_string(index_of_pre_node) + ")_Layer_contain");
	Eigen::Vector3f current_orientation = vector_after.cast<float>();
	if (open_vis_stair_case == true)
		visualize_layers_stair_case(vis_file, vis_file_contain, current_orientation, judge_continue_additive, id_continue);
}

void HybridManufacturing::Visualize_layer_segments(const std::vector<Layer>& layers) const
{
	vector<vector<cv::Point3d>> lines;
	for (const auto& layer : layers) {
		for (const auto& segment : layer.segments) {
			lines.emplace_back();
			auto& segment_line = lines.back();
			segment_line.emplace_back(segment.vertices[0].x, segment.vertices[0].y, layer.z);
			segment_line.emplace_back(segment.vertices[1].x, segment.vertices[1].y, layer.z);
		}
	}

	Visual vis;
	vis.generateModelForRendering_9(lines, file_name + "_katana_layers.obj");
}

void HybridManufacturing::Visualize_layer_polylines(const std::vector<Polylines>& layer_polylines) const
{
	vector<vector<cv::Point3d>> lines;
	for (const auto& layer : layer_polylines) {
		for (const auto& polyline : layer) {
			if (polyline.size() < 2) {
				continue;
			}
			lines.emplace_back();
			auto& polyline_line = lines.back();
			polyline_line.reserve(polyline.size());
			for (const auto& point : polyline) {
				polyline_line.emplace_back(point.x(), point.y(), point.z());
			}
		}
	}

	Visual vis;
	vis.generateModelForRendering_9(lines, file_name + "_cgal_polylines.obj");
}

Layer_Graph HybridManufacturing::BuildAdditiveLayerGraphWithSurfaceMesh(
	const Eigen::Vector3d& vector_after,
	int height_of_beam_search,
	int continue_node_id,
	const nozzle& the_nozzle,
	double& slicing_time,
	double& graph_time) const
{
	clock_t start_time_6 = clock();
	std::vector<SurfaceMeshSliceData> slices;
	BuildSurfaceMeshSlices(slices);

	clock_t end_time_6 = clock();
	slicing_time += double(end_time_6 - start_time_6) / CLOCKS_PER_SEC;

	std::vector<Polylines> contours;
	for (const auto& slice : slices) {
		contours.push_back(slice.contour_points);
		for (auto& contour : contours.back()) {
			contour.push_back(contour.front());
		}
	}
	//Visualize_layer_polylines(contours);
	std::vector<std::vector<std::vector<Eigen::Vector3d>> > surface_slice_points;
	std::vector<std::vector<std::vector<Eigen::Vector3d>> > surface_slice_points_contain;
	surface_slice_points.reserve(slices.size());
	surface_slice_points_contain.reserve(slices.size());
	for (const auto& slice : slices) {
		surface_slice_points.emplace_back();
		surface_slice_points_contain.emplace_back();
		for (const auto& contour : slice.contour_points) {
			std::vector<Eigen::Vector3d> contour_points;
			contour_points.reserve(contour.size());
			for (const auto& point : contour) {
				contour_points.push_back(Eigen::Vector3d(point.x(), point.y(), point.z()));
			}
			surface_slice_points.back().push_back(std::move(contour_points));
			surface_slice_points_contain.back().emplace_back();
		}
	}

	Data surface_data;
	surface_data.ReadData(surface_slice_points, surface_slice_points_contain);
	Layer_Graph layer_graph(surface_data);
	clock_t start_time_4 = clock();
	layer_graph.BuildLayerGraphFromSurfaceMeshSlices(slices, the_nozzle);
	clock_t end_time_4 = clock();
	graph_time += double(end_time_4 - start_time_4) / CLOCKS_PER_SEC;

	return layer_graph;
}

std::vector<std::vector<int>> HybridManufacturing::EvaluateMergedPatchToolCollision(
	const Slicer_2& merged_patch,
	const std::vector<int>& merged_face_source_patch_id,
	cutter cutting_tool,
	bool is_local = false) const
{
	const int face_count = static_cast<int>(merged_patch.triangles.size());
	const int ori_count = static_cast<int>(sampling_subtractive.sample_points.size());

	// face x orientation: -1 表示无碰撞；>=0 表示碰撞面的最大 patch_index
	std::vector<std::vector<int>> min_collision_patch_matrix(face_count, std::vector<int>(ori_count, -1));
	if (face_count == 0 || ori_count == 0) {
		return min_collision_patch_matrix;
	}

	// 每个三角面的采样点（内心）
	std::vector<vasco::core::Vec3> sample_points(face_count);
	for (int i = 0; i < face_count; ++i) {
		const auto& tri = merged_patch.triangles[i];
		const auto& v1 = merged_patch.positions[tri[0]];
		const auto& v2 = merged_patch.positions[tri[1]];
		const auto& v3 = merged_patch.positions[tri[2]];

		double a = distance3d(v1, v2);
		double b = distance3d(v1, v3);
		double c = distance3d(v2, v3);

		vasco::core::Vec3 p{};
		if (a + b + c == 0.0) {
			p = v1;
		}
		else {
			for (int k = 0; k < 3; ++k) {
				p[k] = (a * v1[k] + b * v2[k] + c * v3[k]) / (a + b + c);
			}
		}
		sample_points[i] = p;
	}
	// 每个三角面的采样点（内心）
	Slicer_2 global_slicer;
	std::vector<vasco::core::Vec3> sample_points_global;
	int global_face_count;
	if (is_local) {
		global_slicer.load(file_name + ".obj");
		sample_points_global.resize(global_slicer.triangles.size());
		global_face_count = static_cast<int>(global_slicer.triangles.size());
		for (int i = 0; i < global_face_count; ++i) {
			const auto& tri = global_slicer.triangles[i];
			const auto& v1 = global_slicer.positions[tri[0]];
			const auto& v2 = global_slicer.positions[tri[1]];
			const auto& v3 = global_slicer.positions[tri[2]];

			double a = distance3d(v1, v2);
			double b = distance3d(v1, v3);
			double c = distance3d(v2, v3);

			vasco::core::Vec3 p{};
			if (a + b + c == 0.0) {
				p = v1;
			}
			else {
				for (int k = 0; k < 3; ++k) {
					p[k] = (a * v1[k] + b * v2[k] + c * v3[k]) / (a + b + c);
				}
			}
			sample_points_global[i] = p;
		}
	}


	PrepareToolForCollision(cutting_tool);

	for (int ori = 0; ori < ori_count; ++ori) {
		std::vector<std::vector<Eigen::Vector3d>> temp_faces(face_count, std::vector<Eigen::Vector3d>(3));
		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				int vid = merged_patch.triangles[i][k];
				temp_faces[i][k](0, 0) = merged_patch.positions[vid][0];
				temp_faces[i][k](1, 0) = merged_patch.positions[vid][1];
				temp_faces[i][k](2, 0) = merged_patch.positions[vid][2];
			}
		}

		std::vector<Eigen::Vector3d> temp_samples(face_count);
		for (int i = 0; i < face_count; ++i) {
			temp_samples[i](0, 0) = sample_points[i][0];
			temp_samples[i](1, 0) = sample_points[i][1];
			temp_samples[i](2, 0) = sample_points[i][2];
		}

		std::vector<std::vector<Eigen::Vector3d>> temp_faces_global(global_face_count, std::vector<Eigen::Vector3d>(3));
		std::vector<Eigen::Vector3d> temp_samples_global(global_face_count);

		if (is_local) {
			for (int i = 0; i < global_face_count; ++i) {
				for (int k = 0; k < 3; ++k) {
					int vid = global_slicer.triangles[i][k];
					temp_faces_global[i][k](0, 0) = global_slicer.positions[vid][0];
					temp_faces_global[i][k](1, 0) = global_slicer.positions[vid][1];
					temp_faces_global[i][k](2, 0) = global_slicer.positions[vid][2];
				}
			}
			for (int i = 0; i < global_face_count; ++i) {
				temp_samples_global[i](0, 0) = sample_points_global[i][0];
				temp_samples_global[i](1, 0) = sample_points_global[i][1];
				temp_samples_global[i](2, 0) = sample_points_global[i][2];
			}
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();

		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();

		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				temp_faces[i][k] = rotMatrix.inverse() * temp_faces[i][k];
			}
			temp_samples[i] = rotMatrix.inverse() * temp_samples[i];
		}

		if (is_local) {
			for (int i = 0; i < global_face_count; ++i) {
				for (int k = 0; k < 3; ++k) {
					temp_faces_global[i][k] = rotMatrix.inverse() * temp_faces_global[i][k];
				}
				temp_samples_global[i] = rotMatrix.inverse() * temp_samples_global[i];
			}
		}

		std::vector<double> max_z_of_faces(face_count, MIN_D);
		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				max_z_of_faces[i] = std::max(max_z_of_faces[i], temp_faces[i][k](2, 0));
			}
		}

		std::vector<double> max_z_of_faces_global(global_face_count, MIN_D);
		if (is_local) {
			for (int i = 0; i < global_face_count; ++i) {
				for (int k = 0; k < 3; ++k) {
					max_z_of_faces_global[i] = std::max(max_z_of_faces_global[i], temp_faces_global[i][k](2, 0));
				}
			}
		}

		std::vector<Eigen::Vector3d> normals(face_count);
		for (int i = 0; i < face_count; ++i) {
			const Eigen::Vector3d v1(temp_faces[i][0](0, 0), temp_faces[i][0](1, 0), temp_faces[i][0](2, 0));
			const Eigen::Vector3d v2(temp_faces[i][1](0, 0), temp_faces[i][1](1, 0), temp_faces[i][1](2, 0));
			const Eigen::Vector3d v3(temp_faces[i][2](0, 0), temp_faces[i][2](1, 0), temp_faces[i][2](2, 0));

			Eigen::Vector3d n = (v2 - v1).cross(v3 - v1);
			n.normalize();
			normals[i] = n;
		}

		std::vector<Eigen::Vector3d> normals_global(global_face_count);
		if (is_local) {
			for (int i = 0; i < global_face_count; ++i) {
				const Eigen::Vector3d v1(temp_faces_global[i][0](0, 0), temp_faces_global[i][0](1, 0), temp_faces_global[i][0](2, 0));
				const Eigen::Vector3d v2(temp_faces_global[i][1](0, 0), temp_faces_global[i][1](1, 0), temp_faces_global[i][1](2, 0));
				const Eigen::Vector3d v3(temp_faces_global[i][2](0, 0), temp_faces_global[i][2](1, 0), temp_faces_global[i][2](2, 0));
				Eigen::Vector3d n = (v2 - v1).cross(v3 - v1);
				n.normalize();
				normals_global[i] = n;
			}
		}

		// 计算：刀尖在第 i 面采样点、方向 ori 时，碰撞到的最大 patch_index
		for (int i = 0; i < face_count; ++i) {
			Eigen::Vector3d center_point;
			center_point.x() = temp_samples[i](0, 0) + cutting_tool.cylinder_r * normals[i].x();
			center_point.y() = temp_samples[i](1, 0) + cutting_tool.cylinder_r * normals[i].y();
			center_point.z() = temp_samples[i](2, 0) + cutting_tool.cylinder_r * normals[i].z();

			int max_patch_idx = std::numeric_limits<int>::min();
			bool has_collision = false;

			for (int ii = 0; ii < face_count; ++ii) {
				if (ii == i) {
					continue;
				}
				if (CheckToolCollisionWithCell(center_point, temp_faces[ii], max_z_of_faces[ii], cutting_tool, 30.0, 3.0)) {
					has_collision = true;
					if (ii < static_cast<int>(merged_face_source_patch_id.size())) {
						max_patch_idx = std::max(max_patch_idx, merged_face_source_patch_id[ii]);
					}
				}
			}

			// 无碰撞记为 -1；有碰撞则赋值为patch index + 1
			if (!has_collision) {
				min_collision_patch_matrix[i][ori] = -1;
				if (is_local) {
					for (int ii = 0; ii < global_face_count; ++ii) {
						if (CheckToolCollisionWithCell(center_point, temp_faces_global[ii], max_z_of_faces_global[ii], cutting_tool, 30.0, 3.0)) {
							has_collision = true;
							break;
						}
					}
					if (has_collision) {
						min_collision_patch_matrix[i][ori] = 100;
					}
				}
			}
			else {
				min_collision_patch_matrix[i][ori] = max_patch_idx + 1;
			}
		}
	}

	//std::vector<int> min_collision_patch_per_face(face_count, -1);
	//for (int i = 0; i < face_count; ++i) {
	//	int row_min = std::numeric_limits<int>::max();
	//	for (int ori = 0; ori < ori_count; ++ori) {
	//		row_min = std::min(row_min, min_collision_patch_matrix[i][ori]);
	//	}
	//	if (row_min == std::numeric_limits<int>::max()) {
	//		row_min = -1;
	//	}
	//	min_collision_patch_per_face[i] = row_min;
	//}

	//// 可选输出
	//std::cout << "[Info] min_collision_patch_per_face size="
	//	<< min_collision_patch_per_face.size() << std::endl;

	//// 根据 min_collision_patch_per_face 输出带颜色的 OBJ（每个面一个颜色）
	//ExportMergedPatchFaceColorOBJ(
	//	merged_patch,
	//	min_collision_patch_per_face,
	//	".\\vis\\merged_patch_face_color_by_min_collision_patch.obj");

	return min_collision_patch_matrix;
}

void HybridManufacturing::ExportMergedPatchFaceColorOBJ(
	const Slicer_2& merged_patch,
	const std::vector<int>& max_collision_patch_per_face,
	const std::string& color_obj_file) const
{
	std::ofstream ofs(color_obj_file);
	if (!ofs.is_open()) {
		std::cout << "[Warn] cannot open file for writing: " << color_obj_file << std::endl;
		return;
	}

	int max_nonneg_label = -1;
	for (int v : max_collision_patch_per_face) {
		if (v >= 0) {
			max_nonneg_label = std::max(max_nonneg_label, v);
		}
	}
	std::cout << "[EvaluateMergedPatchToolCollision] max_nonneg_label = "
		<< max_nonneg_label << std::endl;

	auto color_from_label = [max_nonneg_label](int label) -> std::array<double, 3> {
		// -1: 无碰撞，固定颜色（浅灰）
		if (label < 0) {
			return { 0.80, 0.80, 0.80 };
		}

		// 非负: 越大越深（蓝色系）
		double t = 0.0;
		if (max_nonneg_label > 0) {
			t = static_cast<double>(label) / static_cast<double>(max_nonneg_label);
		}
		t = std::max(0.0, std::min(1.0, t));

		const std::array<double, 3> light = { 0.78, 0.88, 1.00 };
		const std::array<double, 3> dark = { 0.05, 0.20, 0.60 };

		std::array<double, 3> c;
		for (int k = 0; k < 3; ++k) {
			c[k] = light[k] * (1.0 - t) + dark[k] * t;
		}
		return c;
		};

	// 每个面写3个独立顶点，保证“面颜色”生效
	int v_count = 0;
	for (int i = 0; i < static_cast<int>(merged_patch.triangles.size()); ++i) {
		const auto& tri = merged_patch.triangles[i];
		const auto color = color_from_label(max_collision_patch_per_face[i]);

		for (int k = 0; k < 3; ++k) {
			const auto& p = merged_patch.positions[tri[k]];
			ofs << "v " << p[0] << " " << p[1] << " " << p[2] << " "
				<< color[0] << " " << color[1] << " " << color[2] << "\n";
		}

		ofs << "f " << (v_count + 1) << " " << (v_count + 2) << " " << (v_count + 3) << "\n";
		v_count += 3;
	}

	std::cout << "[Info] wrote colored OBJ: " << color_obj_file << std::endl;
}

Slicer_2 HybridManufacturing::MergeBlockPatchesWithDedup(
	int max_patch_index,
	std::vector<int>& merged_vertex_source_patch_id,
	std::vector<int>& merged_face_source_patch_id,
	double merge_eps) const
{
	Slicer_2 merged;
	merged_vertex_source_patch_id.clear();
	merged_face_source_patch_id.clear();

	struct QuantizedVertexKey {
		long long x;
		long long y;
		long long z;
		bool operator<(const QuantizedVertexKey& rhs) const {
			if (x != rhs.x) return x < rhs.x;
			if (y != rhs.y) return y < rhs.y;
			return z < rhs.z;
		}
	};

	const double inv_eps = 1.0 / merge_eps;
	auto make_vertex_key = [inv_eps](const vasco::Slicer::Vec3& p) {
		return QuantizedVertexKey{
			static_cast<long long>(std::llround(p[0] * inv_eps)),
			static_cast<long long>(std::llround(p[1] * inv_eps)),
			static_cast<long long>(std::llround(p[2] * inv_eps))
		};
		};

	std::map<QuantizedVertexKey, int> vertex_map;
	std::map<std::array<int, 3>, int> face_key_to_merged_face_index;

	for (int patch_index = 1; patch_index <= max_patch_index; ++patch_index) {
		Slicer_2 patch;
		const std::string patch_file = ".\\vis\\block_patch-" + to_string(patch_index) + "_0" + ".obj";

		if (!patch.load(patch_file)) {
			std::cout << "[Warn] cannot load patch file: " << patch_file << std::endl;
			continue;
		}

		std::vector<int> local_to_global(patch.positions.size(), -1);

		for (int i = 0; i < static_cast<int>(patch.positions.size()); ++i) {
			const auto key = make_vertex_key(patch.positions[i]);
			auto it = vertex_map.find(key);

			if (it == vertex_map.end()) {
				const int new_index = static_cast<int>(merged.positions.size());
				merged.positions.push_back(patch.positions[i]);
				vertex_map.insert({ key, new_index });
				local_to_global[i] = new_index;

				// 单归属：首次出现即最小 patch_index
				merged_vertex_source_patch_id.push_back(patch_index);
			}
			else {
				// 已存在：保持原归属（更小 patch_index）
				local_to_global[i] = it->second;
			}
		}

		for (const auto& tri : patch.triangles) {
			int a = local_to_global[tri[0]];
			int b = local_to_global[tri[1]];
			int c = local_to_global[tri[2]];

			if (a == b || b == c || a == c) {
				continue;
			}

			std::array<int, 3> key = { a, b, c };
			std::sort(key.begin(), key.end());

			auto fit = face_key_to_merged_face_index.find(key);
			if (fit == face_key_to_merged_face_index.end()) {
				const int new_face_index = static_cast<int>(merged.triangles.size());
				merged.triangles.push_back({ a, b, c });
				face_key_to_merged_face_index.insert({ key, new_face_index });

				// 单归属：首次出现即最小 patch_index
				merged_face_source_patch_id.push_back(patch_index);
			}
			// 已存在：保持原归属（更小 patch_index）
		}
	}

	return merged;
}


void HybridManufacturing::CutMesh(
	CutLayerVector all_layers,
	CutLayerVector all_layers_contain,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d vector_after,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	vector<int> all_cut_layers_dependency_layer,
	bool& jud_outer_beam_search_terminate,
	vector<TRiangle>& current_remove_triangles,
	Slicer_2& current_slicer,
	bool judge_continue_additive,
	bool flag_is_continue_block,
	int pre_cont_number_of_queue,
	bool& jud_error,
	int id_node,
	int id_continue,
	vector<int> flag_cut_layers_is_hole)
{
	std::cout << "[HybridManufacturing::CutMesh] height_of_beam_search: " << height_of_beam_search
		<< " cont_number_of_queue: " << cont_number_of_queue
		<< " index_of_pre_node: " << index_of_pre_node
		<< " id_node: " << id_node
		<< " id_continue: " << id_continue
		<< std::endl;
	bool using_solid_model = true;

	std::vector<bool> jud_triangle_have_been_added;	//三角形是否被添加进remove_triangles
	vector<vector<int>> cutting_plane_points;	//每个切割平面对应的切割点
	cutting_plane_points.resize(all_cut_layers.size());
	vector<vector<pair<int, int>>> cutting_plane_edges;	//每个切割平面对应的切割边 <顶点index,顶点index>
	cutting_plane_edges.resize(all_cut_layers.size());

	if (all_cut_layers.size() == 0) {
		std::cout << "[HybridManufacturing::CutMesh]No cut layers!" << std::endl;
		return;
	}
	//////sort cut layers//////
	map<int, int> follow_index;	//记录排序前后cut layers的对应关系 <排序后index,排序前index>
	SortCutLayersByHeight(all_cut_layers, all_cut_layers_dependency_layer, flag_cut_layers_is_hole, follow_index);
	///////////////////////////
	//cout << "a" << id_continue;
	Slicer_2 slicer = LoadSlicerForCutMesh(flag_is_continue_block, height_of_beam_search, index_of_pre_node, id_continue);


	if (flag_is_continue_block == true) {
		height_of_beam_search--;
		cont_number_of_queue = pre_cont_number_of_queue;
	}

	//need rotate first//
	Eigen::Vector3d vectorBefore(0, 0, 1);
	RotateSlicerPositions(slicer, vectorBefore, vector_after);
	/////////////////////

	//layer visualization//
	RotateLayersForVisualization(all_layers, all_layers_contain, vector_after, vectorBefore);
	cout << "b" << endl;
	VisualizeCutLayers(all_layers, all_layers_contain, height_of_beam_search, cont_number_of_queue, index_of_pre_node, judge_continue_additive, id_continue, vector_after);
	////////////////////////
	cout << "B" << endl;
	//--------------------cut-----------------------//
	slicer.normal[0] = 0;
	slicer.normal[1] = 0;
	slicer.normal[2] = 1;

	slicer.jud_plane.resize(slicer.triangles.size());
	for (int i = 0; i < slicer.triangles.size(); i++)
		slicer.jud_plane[i] = false;
	clock_t start_time_3, end_time_3;

	for (int i = 0; i < all_cut_layers.size(); i++) {
		slicer.origin[0] = 0;
		slicer.origin[1] = 0;
		slicer.origin[2] = all_cut_layers[i][0].z;
		slicer.cut();
		//cout << "ok";
		cout << "i " << i << " slicer.origin[2] " << slicer.origin[2] << endl;
	}

	std::vector<TRiangle> ori_triangle = slicer.triangles;	//保存切割后的所有三角形
	jud_triangle_have_been_added.resize(slicer.triangles.size());	//为ori_triangle中的每个三角形分配一个标记，表示是否被添加进remove_triangles
	for (int i = 0; i < slicer.triangles.size(); i++)
		jud_triangle_have_been_added[i] = false;

	Slicer_2 all_slicer;	//all_slicer保存切割后的所有三角形和顶点

	all_slicer.positions = slicer.positions;
	all_slicer.triangles = ori_triangle;


	all_slicer.save(file_name + "-after_cut-" + std::to_string(height_of_beam_search) + "EEEERQ.obj");

	//--------------------save candidate_triangles-----------------------//
	int current_index = 0;	//
	std::vector<VEctor> min_z_point;	//每个候选三角形中z值最小的顶点坐标
	std::vector<double> min_z_triangle;	//每个候选三角形中z值最小的顶点的z值
	vector<int> index_of_min_point_in_triangle;	//每个候选三角形中z值最小的顶点在三角形中的索引
	std::vector<TRiangle> candidate_triangles;	//候选三角形集合
	std::vector<int> id_candidate_triangles;	//候选三角形在all_slicer.triangles中的索引
	std::vector<int> id_triangles;	//候选三角形在slicer.triangles中的索引
	std::vector<TRiangle> remove_triangles;	//最终被移除的三角形集合
	std::vector<int> id_remove_triangles;	//最终被移除的三角形在all_slicer.triangles中的索引
	candidate_triangles.clear();
	id_triangles.clear();
	double min_z_all_cut_layers = 9999999;	//所有切割层中z值最小的值
	for (int t = 0; t < all_cut_layers.size(); t++)
		min_z_all_cut_layers = min(min_z_all_cut_layers, all_cut_layers[t][0].z);
	start_time_3 = clock();
	for (int i = 0; i < all_slicer.triangles.size(); i++) {	//遍历所有切割后的三角形
		double temp_min_z_triangle = 9999999;	//当前三角形中z值最小的顶点的z值
		for (int k = 0; k < 3; k++) {
			temp_min_z_triangle = min(all_slicer.positions[all_slicer.triangles[i][k]][2], temp_min_z_triangle);
		}
		if (temp_min_z_triangle + 0.001 >= min_z_all_cut_layers) {	//如果当前三角形中z值最小的顶点的z值大于等于所有切割层中z值最小的值，则将该三角形加入候选三角形集合
			candidate_triangles.push_back(all_slicer.triangles[i]);	//加入候选三角形集合
			id_candidate_triangles.push_back(i);	//记录该候选三角形在all_slicer.triangles中的索引
			for (int k = i; k < slicer.triangles.size(); k++)
				if (slicer.triangles[k] == all_slicer.triangles[i]) {	//找到该候选三角形在slicer.triangles中的索引
					id_triangles.push_back(k);
					break;
				}
		}
	}
	end_time_3 = clock();
	for (int i = 0; i < candidate_triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;	//当前候选三角形中z值最小的顶点的z值
		VEctor temp_min_z_point = all_slicer.positions[candidate_triangles[i][0]];
		int temp_index_of_min_point_in_triangle;
		for (int k = 0; k < 3; k++) {
			if (all_slicer.positions[candidate_triangles[i][k]][2] < temp_min_z_triangle) {	//找到当前候选三角形中z值最小的顶点
				temp_min_z_triangle = all_slicer.positions[candidate_triangles[i][k]][2];	//更新z值最小的顶点的z值
				temp_min_z_point = all_slicer.positions[candidate_triangles[i][k]];	//更新z值最小的顶点的坐标
				temp_index_of_min_point_in_triangle = k;	//更新z值最小的顶点在三角形中的索引
			}
		}
		min_z_triangle.push_back(temp_min_z_triangle);	//将当前候选三角形中z值最小的顶点的z值加入min_z_triangle
		min_z_point.push_back(temp_min_z_point);	//将当前候选三角形中z值最小的顶点的坐标加入min_z_point
		index_of_min_point_in_triangle.push_back(temp_index_of_min_point_in_triangle);	//将当前候选三角形中z值最小的顶点在三角形中的索引加入index_of_min_point_in_triangle
	}

	cout << "()()()(" << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << endl;
	clock_t start_time_2, end_time_2;
	start_time_2 = clock();
	for (int i = 0; i < candidate_triangles.size(); i++) {	//对候选三角形按照min_z_triangle进行排序
		for (int j = i + 1; j < candidate_triangles.size(); j++) {
			if (min_z_triangle[i] > min_z_triangle[j]) {
				swap(candidate_triangles[i], candidate_triangles[j]);	//交换候选三角形
				swap(id_candidate_triangles[i], id_candidate_triangles[j]);	//交换候选三角形在all_slicer.triangles中的索引
				swap(id_triangles[i], id_triangles[j]);	//交换候选三角形在slicer.triangles中的索引
				swap(min_z_triangle[i], min_z_triangle[j]);	//交换min_z_triangle
				swap(min_z_point[i], min_z_point[j]);	//交换min_z_point
				swap(index_of_min_point_in_triangle[i], index_of_min_point_in_triangle[j]);	//交换index_of_min_point_in_triangle
			}
		}
	}
	end_time_2 = clock();
	cout << "()()()(" << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << endl;
	//slicer.triangles = candidate_triangles;

	//--------------------save OPP_triangles one by one-----------------------//
	vector<int> save_current_index;
	for (int t = 0; t < all_cut_layers.size(); t++) {	//枚举all_cut_layers中的每一层
		current_index = 0;   //该方式也许比较慢
		double boundary_bottom = 999999, boundary_left = 999999, boundary_top = -999999, boundary_right = -999999;	//当前切割层的边界
		//Point_2* points = new Point_2[all_cut_layers[t].size()];
		vector<Point_2> points;
		points.resize(all_cut_layers[t].size());
		for (int i = 0; i < all_cut_layers[t].size(); i++) {	//枚举all_cut_layers[t]中的点，计算当前切割层的二维AABB
			boundary_top = std::max(boundary_top, all_cut_layers[t][i].y);
			boundary_bottom = std::min(boundary_bottom, all_cut_layers[t][i].y);
			boundary_right = std::max(boundary_right, all_cut_layers[t][i].x);
			boundary_left = std::min(boundary_left, all_cut_layers[t][i].x);
			Point_2 temp_point(all_cut_layers[t][i].x, all_cut_layers[t][i].y);
			points[i] = temp_point;
		}
		cv::Point2d current_triangle_point;
		cv::Point2d current_layer_point;
		for (; current_index < candidate_triangles.size(); current_index++) {	//枚举所有候选三角形
			if (abs(min_z_triangle[current_index] - all_cut_layers[t][0].z) > 0.0001) {
				if (min_z_triangle[current_index] > all_cut_layers[t][0].z) {	//如果当前候选三角形中z值最小的顶点的z值大于当前切割层的z值，则跳出循环
					break;
				}
			}
			else {	//如果当前候选三角形中最小z值等于当前切割层的z值,才进行是否加入remove_triangles的判断
				bool jud_is_boundary_point = false;
				//for (int i = 0; i < all_cut_layers[t].size(); i++) {
				int cont_inside_boundary = 0;
				for (int k = 0; k < 3; k++) {	//枚举当前候选三角形的三个顶点，判断是否在当前切割层的AABB内
					if (all_slicer.positions[candidate_triangles[current_index][k]][0] + 0.01 >= boundary_left && all_slicer.positions[candidate_triangles[current_index][k]][0] - 0.01 <= boundary_right
						&& all_slicer.positions[candidate_triangles[current_index][k]][1] + 0.01 >= boundary_bottom && all_slicer.positions[candidate_triangles[current_index][k]][1] - 0.01 <= boundary_top) {
						cont_inside_boundary++;	//判断当前顶点在AABB内，计数加1
					}
				}
				if (cont_inside_boundary >= 2)
					jud_is_boundary_point = true;
				/*if (min_z_point[current_index][0] + 0.01 >= boundary_left && min_z_point[current_index][0] - 0.01 <= boundary_right &&
					min_z_point[current_index][1] + 0.01 >= boundary_bottom && min_z_point[current_index][1] - 0.01 <= boundary_top) {
					jud_is_boundary_point = true;
					break;
				}*/
				//}

				if (jud_is_boundary_point == true) {
					auto layer_begin = points.begin();
					auto layer_end = points.end();
					//if (check_inside_2(Point_2(current_triangle_point.x, current_triangle_point.y), points, points + all_cut_layers[t].size(), K())) {
					current_triangle_point.x = min_z_point[current_index][0];
					current_triangle_point.y = min_z_point[current_index][1];
					/*if (check_inside_2(Point_2(current_triangle_point.x, current_triangle_point.y), layer_begin, layer_end, K())) {*/

					if (1) {
						int cont_cutting_points = 0;
						int temp_left = 0, temp_right = 0;
						for (int k = 0; k < 3; k++) {
							if (abs(all_cut_layers[t][0].z - all_slicer.positions[candidate_triangles[current_index][k]][2]) < 0.0001) {	//判断当前顶点是否在切割平面上（前面已经保证了z值最小的顶点不大于切割平面z值）
								cont_cutting_points++;
								if (cont_cutting_points == 1)
									temp_left = candidate_triangles[current_index][k];
								else if (cont_cutting_points == 2) {
									temp_right = candidate_triangles[current_index][k];
									cutting_plane_edges[t].push_back(make_pair(temp_left, temp_right));
									break;
								}
								//cutting_plane_points[t].push_back(candidate_triangles[current_index][index_of_min_point_in_triangle[current_index]]);
							}
						}
						remove_triangles.push_back(candidate_triangles[current_index]);	//将当前候选三角形加入remove_triangles
						id_remove_triangles.push_back(id_candidate_triangles[current_index]);	//将当前候选三角形在all_slicer.triangles中的索引加入id_remove_triangles
						/*if (height_of_beam_search == 2 && id_triangles[current_index] == 10297)
							cout << "********************************** " << all_cut_layers[t][0].z << endl;*/

						save_current_index.push_back(current_index);	//记录当前候选三角形在candidate_triangles中的索引到save_current_index
						jud_triangle_have_been_added[id_triangles[current_index]] = true;	//jud_triangle_have_been_added标记当前候选三角形已被添加进remove_triangles
					}
					else {
						for (int j = 0; j < all_cut_layers[t].size(); j++) {
							current_layer_point.x = all_cut_layers[t][j].x;	//
							current_layer_point.y = all_cut_layers[t][j].y;
							current_triangle_point.x = min_z_point[current_index][0];
							current_triangle_point.y = min_z_point[current_index][1];
							//if ((jud_is_boundary_point == false && distance_2d(current_layer_point, current_triangle_point) < 0.002) || (jud_is_boundary_point == true)) {   //0.002
							if (distance2d(current_layer_point, current_triangle_point) < 0.1) {  //4.0
								int cont_cutting_points = 0;
								int temp_left = 0, temp_right = 0;
								for (int k = 0; k < 3; k++) {
									if (abs(all_cut_layers[t][0].z - all_slicer.positions[candidate_triangles[current_index][k]][2]) < 0.0001) {
										cont_cutting_points++;
										if (cont_cutting_points == 1)
											temp_left = candidate_triangles[current_index][k];
										else if (cont_cutting_points == 2) {
											temp_right = candidate_triangles[current_index][k];
											cutting_plane_edges[t].push_back(make_pair(temp_left, temp_right));
											break;
										}
										//cutting_plane_points[t].push_back(candidate_triangles[current_index][index_of_min_point_in_triangle[current_index]]);
									}
								}
								remove_triangles.push_back(candidate_triangles[current_index]);
								id_remove_triangles.push_back(id_candidate_triangles[current_index]);
								/*if (height_of_beam_search == 2 && id_triangles[current_index] == 10297)
									cout << "********************************** " << all_cut_layers[t][0].z << endl;*/

								save_current_index.push_back(current_index);
								jud_triangle_have_been_added[id_triangles[current_index]] = true;
								break;
							}
						}
					}
				}
			}

		}
	}
	current_index = 0;

	while (1) {
		bool flag_break = false;
		for (int i = 0; i < remove_triangles.size(); i++)
		{
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					if (jud_triangle_have_been_added[id_triangles[current_index]] == false && candidate_triangles[current_index][j] == remove_triangles[i][k]
						&& min_z_triangle[current_index] >= min_z_triangle[save_current_index[i]]) {	//如果当前候选三角形的某个顶点与remove_triangles[i]的某个顶点相同，且当前候选三角形中z值最小的顶点的z值大于等于remove_triangles[i]中z值最小的顶点的z值
						//if (height_of_beam_search != 2)
						remove_triangles.push_back(candidate_triangles[current_index]);	//将当前候选三角形加入remove_triangles
						id_remove_triangles.push_back(id_candidate_triangles[current_index]);
						save_current_index.push_back(current_index);
						jud_triangle_have_been_added[id_triangles[current_index]] = true;
						flag_break = true;
						break;
					}
				}
				if (flag_break == true)
					break;
			}
			if (flag_break == true)
				break;
		}
		current_index++;
		if (current_index >= candidate_triangles.size())
			break;
	}

	vasco::Slicer temp_slicer_1;
	temp_slicer_1.clear();
	temp_slicer_1.positions = all_slicer.positions;
	for (int i = 0; i < id_remove_triangles.size(); i++) {
		temp_slicer_1.triangles.push_back(all_slicer.triangles[id_remove_triangles[i]]);
	}
	//cout << "remove faces " << temp_slicer_1.triangles.size() << endl;

	//RotateSlicerPositions(temp_slicer_1, vector_after, vectorBefore);
	temp_slicer_1.save(file_name + "-all2-" + std::to_string(height_of_beam_search) + "QWERQ.obj");
	
	std::vector<std::vector<Eigen::Vector3d>> all_cut_layers_v3d;
	for (int i = 0; i < all_cut_layers.size(); i++) {
		std::vector<Eigen::Vector3d> temp_layer;
		for (int j = 0; j < all_cut_layers[i].size(); j++) {
			Eigen::Vector3d temp_point(all_cut_layers[i][j].x, all_cut_layers[i][j].y, all_cut_layers[i][0].z);
			temp_layer.push_back(temp_point);
		}
		all_cut_layers_v3d.push_back(temp_layer);
	}
	vasco::io::outputContourToObj(all_cut_layers_v3d, Eigen::Vector3d(0.5, 0.5, 0.0),  file_name + "-all_cut_layers-" + std::to_string(height_of_beam_search) + "QWERQ.obj");


	//add remaining face, set a distance threshold Dis, only a face less Dis from other cut layers and no other dependency layer exist, the face is consider to remaining face
	double Dis = dh * 2;
	for (int i = 0; i < jud_triangle_have_been_added.size(); i++) {
		if (jud_triangle_have_been_added[i] == false) {
			int id_layer;
			double min_dis = 99999999;
			for (int t = 0; t < all_cut_layers.size(); t++) {
				if (all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] <= 2 * dh && all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] >= -0.001) {	//如果当前切割层的z值与当前三角形中z值最小的顶点的z值之差在[-0.001,2*dh]范围内
					for (int j = 0; j < all_cut_layers[t].size(); j++) {
						cv::Point3d current_triangle_point(slicer.positions[slicer.triangles[i][0]][0], slicer.positions[slicer.triangles[i][0]][1], slicer.positions[slicer.triangles[i][0]][2]);
						cv::Point3d current_layer_point(all_cut_layers[t][j].x, all_cut_layers[t][j].y, all_cut_layers[t][0].z);

						double distance = distance3d(current_triangle_point, current_layer_point);
						if (distance < min_dis) {
							min_dis = distance;
							id_layer = t;
						}
					}
				}
			}
			if (min_dis < Dis && all_cut_layers_dependency_layer[id_layer] == 0) {
				cout << "还真的能到这个地方啊 - 残余面片 min_dis " << min_dis << "Dis " << Dis << std::endl;
				remove_triangles.push_back(slicer.triangles[i]);
				id_remove_triangles.push_back(i);
				//jud_triangle_have_been_added[i] = true;
			}

		}
	}




	/////////////////////////删除残余面片////////////////////////////
	//建立面片邻接关系
	vector<vector<int>> adjacent_faces(slicer.triangles.size());
	//建立不在remove_triangles中的面片之间的邻接关系
	for (int i = 0; i < slicer.triangles.size(); i++) {
		if (jud_triangle_have_been_added[i] == false) {
			for (int j = i + 1; j < slicer.triangles.size(); j++) {
				if (jud_triangle_have_been_added[j] == false) {
					int cont_same_point = 0;
					for (int k = 0; k < 3; k++) {
						for (int l = 0; l < 3; l++) {
							if (slicer.triangles[i][k] == slicer.triangles[j][l]) {
								cont_same_point++;
								break;
							}
						}
					}
					if (cont_same_point >= 2) {	//如果两个面片有两个以上的公共顶点，则认为它们是邻接的
						adjacent_faces[i].push_back(j);
						adjacent_faces[j].push_back(i);
					}
				}
			}
		}
	}
	//分区
	vector<bool> visited(slicer.triangles.size(), false);
	vector<vector<int>> save_faces(0);
	for (int i = 0; i < slicer.triangles.size(); ++i)
	{
		if (jud_triangle_have_been_added[i] || visited[i]) {
			continue; // 已标记或已访问，跳过
		}

		visited[i] = true;

		// 当前连通分量的索引集合
		std::vector<int> component;
		component.push_back(i);

		// 标准队列进行 BFS
		std::queue<int> q;
		q.push(i);
		while (!q.empty())
		{
			const int u = q.front();
			q.pop();

			// 遍历 u 的邻接面
			const auto& adj = adjacent_faces[u];
			for (int v : adj)
			{
				if (!visited[v] && !jud_triangle_have_been_added[v])
				{
					visited[v] = true;
					q.push(v);
					component.push_back(v);
				}
			}
		}

		// 将本连通分量保存（对应之前的 save_faces.push_back(...)）
		save_faces.push_back(std::move(component));
	}

	// 移除过小的残余面片分量
	for (const auto& comp : save_faces)
	{
		cout << "&&*&& " << comp.size() << endl;
		if (comp.size() < 80)
		{
			for (int idx : comp)
			{
				remove_triangles.push_back(slicer.triangles[idx]);
				id_remove_triangles.push_back(idx);
				jud_triangle_have_been_added[idx] = true;
			}
		}
	}
	//////////////////////////////////////////////////////////////////



	cout << "c" << endl;
	//current_remove_triangles = remove_triangles;

	current_remove_triangles = remove_triangles; // current_remove_triangles先不进行筛选，直接等于remove_triangles，后面再剔除非表面的面片
	cout << "ccc" << endl;

	//将slicer旋转回原始位置
	RotateSlicerPositions(slicer, vector_after, vectorBefore);
	current_slicer = slicer;

	//从all_slicer中移除remove_triangles
	//cout << current_remove_triangles.size() << endl;
	//if (height_of_beam_search != 2 || id_continue != 1)
	for (int i = 0; i < remove_triangles.size(); i++) {
		for (int j = 0; j < all_slicer.triangles.size(); ) {
			if (remove_triangles[i][0] == all_slicer.triangles[j][0] && remove_triangles[i][1] == all_slicer.triangles[j][1] && remove_triangles[i][2] == all_slicer.triangles[j][2]) {
				all_slicer.triangles.erase(all_slicer.triangles.begin() + j);
				break;
			}
			j++;
		}
	}	//潜在问题？：如果顶点顺序不同（例如顺时针/逆时针），将无法删除；可改为排序比较或集合比较。
	//这个潜在问题应该不会出现，因为计算remove_triangles时，是直接从all_slicer.triangles中取出的三角形，三角形顶点顺序没有改变。


	size_t face_cnt = Normals.rows();
	std::cout << "face_cnt " << face_cnt << std::endl;
	current_remove_triangles = FilterSurfaceRemoveTriangles(slicer, remove_triangles);

	//add cutting plane triangles
	/*all_slicer.triangles.insert(all_slicer.triangles.begin(),cutting_plane_points)
	cutting_plane_points*/

	//erase some cutting plane
	for (int i = 0; i < all_cut_layers.size(); i++) {
		if (all_cut_layers_dependency_layer[i] == 0) {
			all_cut_layers.erase(all_cut_layers.begin() + i);
			all_cut_layers_dependency_layer.erase(all_cut_layers_dependency_layer.begin() + i);
			cutting_plane_edges.erase(cutting_plane_edges.begin() + i);
			cutting_plane_points.erase(cutting_plane_points.begin() + i);
			i--;
		}
	}
	cout << "cccc" << endl;
	//sort cutting_plane_points by adjacency relation
	vector<vector<int>> real_cutting_plane_triangles;	//为每个切割层存储排序后的切割平面顶点索引
	real_cutting_plane_triangles.resize(all_cut_layers.size());	//real_cutting_plane_triangles大小与all_cut_layers相同
	for (int i = 0; i < all_cut_layers.size(); i++) {
		cv::Point2d current_triangle_point;
		cv::Point2d current_layer_point;
		current_layer_point.x = all_cut_layers[i][0].x;
		current_layer_point.y = all_cut_layers[i][0].y;
		double min_dis = MAX_D;
		int index_start_point_id;
		for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
			current_triangle_point.x = all_slicer.positions[cutting_plane_edges[i][j].second][0];
			current_triangle_point.y = all_slicer.positions[cutting_plane_edges[i][j].second][1];
			if (distance2d(current_layer_point, current_triangle_point) < min_dis)
			{
				index_start_point_id = j;
				min_dis = distance2d(current_layer_point, current_triangle_point);
			}
		}
		int index_start_point = cutting_plane_edges[i][index_start_point_id].second;
		int index_of_last_edge = index_start_point_id;
		real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][index_start_point_id].first);
		//cout << "kkkk" << cutting_plane_edges[i].size() << endl;
		int jud_select_id = -1;
		int cont_segment = 0;
		while (index_start_point != cutting_plane_edges[i][index_start_point_id].first) {
			bool jud_select = false;
			cont_segment++;
			for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
				if (j != index_of_last_edge && cutting_plane_edges[i][j].first == index_start_point) {
					real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].first);
					index_start_point = cutting_plane_edges[i][j].second;
					index_of_last_edge = j;
					jud_select = true;
					jud_select_id = 0;
					break;
				}
				else if (j != index_of_last_edge && cutting_plane_edges[i][j].second == index_start_point) {
					real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].second);
					index_start_point = cutting_plane_edges[i][j].first;
					index_of_last_edge = j;
					jud_select = true;
					jud_select_id = 1;
					break;
				}
			}
			if (cont_segment > 50000) {
				jud_select = false;
				break;
				//real_cutting_plane_triangles[i].clear();
			}
			if (cont_segment > 1000000) {
				jud_select = false;
				//real_cutting_plane_triangles[i].clear();
			}
			if (jud_select == false) {
				//cout << "no:" << i << "  ";
				jud_error = true;
				/*for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
					if (j != index_of_last_edge) {
						cv::Point2d current_point_2, current_point_3, current_point_4;
						current_point_3.x = all_slicer.positions[cutting_plane_edges[i][j].first][0];
						current_point_3.y = all_slicer.positions[cutting_plane_edges[i][j].first][1];
						current_point_4.x = all_slicer.positions[cutting_plane_edges[i][j].second][0];
						current_point_4.y = all_slicer.positions[cutting_plane_edges[i][j].second][1];
						if (jud_select_id == 0) {
							current_point_2.x = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].second][0];
							current_point_2.y = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].second][1];
						}
						else if (jud_select_id == 1) {
							current_point_2.x = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].first][0];
							current_point_2.y = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].first][1];
						}
						if (distance_2d(current_point_2, current_point_3) < 0.00001) {
							real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].first);
							index_start_point = cutting_plane_edges[i][j].second;
							index_of_last_edge = j;
							jud_select = true;
							jud_select_id = 0;
							break;
						}
						else if (distance_2d(current_point_2, current_point_4) < 0.00001) {
							real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].second);
							index_start_point = cutting_plane_edges[i][j].first;
							index_of_last_edge = j;
							jud_select = true;
							jud_select_id = 1;
							break;
						}
					}
				}*/
				break;
			}
		}
	}
	/*if (jud_error[id_node] == true)
		return;*/
		/*cout << "mm" << endl;*/
		/*if (height_of_beam_search == 6)
			return;*/
	for (int i = 0; i < real_cutting_plane_triangles.size(); i++) {
		if (real_cutting_plane_triangles[i].size() == 1) {
			real_cutting_plane_triangles.erase(real_cutting_plane_triangles.begin() + i);
			i--;
		}
	}

	///////////////////////////terminate//////////////////////////////
	if (all_slicer.triangles.size() < terminate_threshold_of_number_of_faces) {
		jud_outer_beam_search_terminate = true;
	}
	//////////////////////////////////////////////////////////////////

	vector<int> id_contact_faces;
	//if (height_of_beam_search != 2 || id_continue != 1)
	if (using_solid_model == true) {
		Anticlockwise(real_cutting_plane_triangles, all_slicer);
		std::map<QuantizedPointKey, int> point_index_map;
		for (const auto& tri : all_slicer.triangles) {
			for (int k = 0; k < 3; ++k) {
				point_index_map.emplace(MakeQuantizedKey(Point_3(all_slicer.positions[tri[k]][0], all_slicer.positions[tri[k]][1], all_slicer.positions[tri[k]][2])), tri[k]);
			}
		}
		for (int i = 0; i < real_cutting_plane_triangles.size(); i++) {
			if (flag_cut_layers_is_hole[i] != -1)
				continue;
			std::vector<int> hole_ring;
			std::vector<std::vector<int>> hole_rings;
			for (int m = 0; m < real_cutting_plane_triangles.size(); m++) {
				if (flag_cut_layers_is_hole[m] != -1 && follow_index[flag_cut_layers_is_hole[m]] == i) {
					hole_ring = real_cutting_plane_triangles[m];
					hole_rings.push_back(hole_ring);
					break;
				}
			}

			contact_face_triangulation_mode = ContactFaceTriangulationMode::Earcut;
			//contact_face_triangulation_mode = ContactFaceTriangulationMode::CGALRemesh;
			if (contact_face_triangulation_mode == ContactFaceTriangulationMode::CGALRemesh) {
				auto generated = TriangulateContactFacesWithCDT(
					real_cutting_plane_triangles[i],
					hole_rings,
					all_slicer,
					id_contact_faces,
					point_index_map);
				continue;
			}

			using Coord = double;
			using NN = uint32_t;
			using Point = std::array<Coord, 2>;
			std::vector<std::vector<Point>> polygon;
			polygon.clear();
			std::vector<Point> temp_vec;
			polygon.push_back(temp_vec);
			map<int, int> map_index_faces;
			map_index_faces.clear();
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				polygon[0].push_back({ all_slicer.positions[real_cutting_plane_triangles[i][j]][0], all_slicer.positions[real_cutting_plane_triangles[i][j]][1] });
				map_index_faces.insert({ j,real_cutting_plane_triangles[i][j] });
			}
			for (int m = 0; m < real_cutting_plane_triangles.size(); m++) {
				if (flag_cut_layers_is_hole[m] != -1 && follow_index[flag_cut_layers_is_hole[m]] == i) {
					for (int j = 0; j < real_cutting_plane_triangles[m].size(); j++) {
						polygon.push_back(temp_vec);
						polygon[1].push_back({ all_slicer.positions[real_cutting_plane_triangles[m][j]][0], all_slicer.positions[real_cutting_plane_triangles[m][j]][1] });
						map_index_faces.insert({ polygon[0].size() + j,real_cutting_plane_triangles[m][j] });
					}
					break;
				}
			}

			std::vector<std::vector<Eigen::Vector3d>> polygons_v3d;
			for (int j = 0; j < polygon.size(); j++) {
				std::vector<Eigen::Vector3d> temp_polygon_v3d;
				for (int k = 0; k < polygon[j].size(); k++) {
					Eigen::Vector3d temp_point_v3d(polygon[j][k][0], polygon[j][k][1], all_slicer.positions[real_cutting_plane_triangles[i][0]][2]);
					temp_polygon_v3d.push_back(temp_point_v3d);
				}
				polygons_v3d.push_back(temp_polygon_v3d);
			}
			vasco::io::outputContourToObj(polygons_v3d, Eigen::Vector3d(0.5, 0.5, 0.0), file_name + "-contact_face_polygon-" + std::to_string(height_of_beam_search) + "_" + std::to_string(i) + ".obj");

			std::vector<NN> indices = mapbox::earcut<NN>(polygon);
			for (int j = 0; j < indices.size();) {
				TRiangle the_new_cutting_plane_triangle;
				the_new_cutting_plane_triangle[0] = map_index_faces[indices[j]]; j++;
				the_new_cutting_plane_triangle[1] = map_index_faces[indices[j]]; j++;
				the_new_cutting_plane_triangle[2] = map_index_faces[indices[j]]; j++;
				all_slicer.triangles.insert(all_slicer.triangles.end(), the_new_cutting_plane_triangle);
				id_contact_faces.push_back(all_slicer.triangles.size() - 1);
			}

		}
		//}
	}




	/*if (height_of_beam_search == 1) {
		for (int i = 0; i < real_cutting_plane_triangles.size(); i++) {
			VEctor new_point;
			double temp_x = 0, temp_y = 0;
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				temp_x += all_slicer.positions[real_cutting_plane_triangles[i][j]][0];
				temp_y += all_slicer.positions[real_cutting_plane_triangles[i][j]][1];
			}
			new_point[0] = temp_x / real_cutting_plane_triangles[i].size();
			new_point[1] = temp_y / real_cutting_plane_triangles[i].size();
			new_point[2] = all_slicer.positions[real_cutting_plane_triangles[i][0]][2];
			all_slicer.positions.insert(all_slicer.positions.end(), new_point);
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				TRiangle the_new_cutting_plane_triangle;
				the_new_cutting_plane_triangle[0] = all_slicer.positions.size() - 1;
				the_new_cutting_plane_triangle[1] = real_cutting_plane_triangles[i][j];
				the_new_cutting_plane_triangle[2] = real_cutting_plane_triangles[i][(j + 1) % real_cutting_plane_triangles[i].size()];
				all_slicer.triangles.insert(all_slicer.triangles.end(), the_new_cutting_plane_triangle);
			}
		}
	}*/

	Slicer_2 all_slicer_2 = all_slicer;
	all_slicer_2.triangles = remove_triangles;
	
	RotateSlicerPositions(all_slicer_2, vector_after, vectorBefore);


	if (judge_continue_additive == false)
		all_slicer_2.save(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_current" + ".obj");
	else
		all_slicer_2.save(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_current" + "_subblock.obj");

	//std::vector<double> cut_plane_z_values;
	//cut_plane_z_values.reserve(all_cut_layers.size());
	//for (const auto& cut_layer : all_cut_layers) {
	//	if (!cut_layer.empty()) {
	//		cut_plane_z_values.push_back(cut_layer[0].z);
	//	}
	//}
	//RemeshCutPatchBeforeSaving(
	//	all_slicer,
	//	id_contact_faces,
	//	cut_plane_z_values,
	//	std::max(dh, 5.0));


	all_slicer.save(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_before_rotate_back.obj");

	RotateSlicerPositions(all_slicer, vector_after, vectorBefore);

	double min_z = MAX_D, max_z = -MAX_D;
	for (int i = 0; i < all_slicer.positions.size(); i++) {
		double zz = all_slicer.positions[i][2];
		if (zz < min_z) {
			min_z = zz;
		}
		if (zz > max_z) {
			max_z = zz;
		}
	}

	if (max_z - min_z < 3.0) {
		jud_outer_beam_search_terminate = true;
	}


	Eigen::MatrixXd temp_V(all_slicer.positions.size(), 3);
	Eigen::MatrixXi temp_F(all_slicer.triangles.size(), 3);
	for (int i = 0; i < all_slicer.positions.size(); i++) {
		temp_V.row(i).x() = all_slicer.positions[i][0];
		temp_V.row(i).y() = all_slicer.positions[i][1];
		temp_V.row(i).z() = all_slicer.positions[i][2];
	}
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		temp_F.row(i).x() = all_slicer.triangles[i][0];
		temp_F.row(i).y() = all_slicer.triangles[i][1];
		temp_F.row(i).z() = all_slicer.triangles[i][2];
	}
	if (judge_continue_additive == false) {
		Katana::Instance().stl.saveStlFromObj(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl", temp_V, temp_F);
		vector<std::array<double, 3>> V3;
		vector<std::array<int, 3>> F3;
		vector<std::array<double, 3>> N3;
		ifstream ifs(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl");
		igl::read_stl_ascii(ifs, V3, F3, N3);
		Eigen::MatrixXd V4(V3.size(), 3);
		Eigen::MatrixXi F4(F3.size(), 3);
		Eigen::MatrixXd N4(N3.size(), 3);
		for (int i = 0; i < V3.size(); i++)
			for (int j = 0; j < 3; j++)
				V4(i, j) = V3[i][j];
		for (int i = 0; i < F3.size(); i++)
			for (int j = 0; j < 3; j++)
				F4(i, j) = F3[i][j];
		igl::writeSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_B.stl", V4, F4, igl::FileEncoding::Binary);

		Geometry tessel;
		tessel.visit(ImportSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_B.stl"));
		tessel.visit(ExportOBJ(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".obj"));
	}
	else {
		Katana::Instance().stl.saveStlFromObj(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock.stl", temp_V, temp_F);
		vector<std::array<double, 3>> V3;
		vector<std::array<int, 3>> F3;
		vector<std::array<double, 3>> N3;
		ifstream ifs(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock.stl");
		igl::read_stl_ascii(ifs, V3, F3, N3);
		Eigen::MatrixXd V4(V3.size(), 3);
		Eigen::MatrixXi F4(F3.size(), 3);
		Eigen::MatrixXd N4(N3.size(), 3);
		for (int i = 0; i < V3.size(); i++)
			for (int j = 0; j < 3; j++)
				V4(i, j) = V3[i][j];
		for (int i = 0; i < F3.size(); i++)
			for (int j = 0; j < 3; j++)
				F4(i, j) = F3[i][j];
		igl::writeSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock_B.stl", V4, F4, igl::FileEncoding::Binary);

		Geometry tessel;
		tessel.visit(ImportSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock_B.stl"));
		tessel.visit(ExportOBJ(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock.obj"));
	}

	string str_contact_faces;
	if (judge_continue_additive == false)
		str_contact_faces = file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_contact_faces.txt";
	else
		str_contact_faces = file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock_contact_faces.txt";

	Visual vis_2;
	cv::Point3d input_ori(vector_after.x(), vector_after.y(), vector_after.z());
	if (judge_continue_additive == false)
		vis_2.generateArrows(input_ori, file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue));
	else
		vis_2.generateArrows(input_ori, file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock");

	ofstream ofile_contact_faces(str_contact_faces);
	ofile_contact_faces << id_contact_faces.size() << endl;
	for (int i = 0; i < id_contact_faces.size(); i++)
		ofile_contact_faces << id_contact_faces[i] << endl;
	cout << "d" << endl;
	/////////////////////////
	return;
}

void HybridManufacturing::subtractive_accessibility_decomposition(
	vector<TRiangle> need_detect_triangle,
	int height_of_beam_search,
	int cont_number_of_queue,
	cutter cutting_tool,
	Slicer_2 current_slicer)
{
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 27;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 23;
	cutting_tool.carriage_height = 33;

	/////show all accessible points in every orientation/////
	Eigen::MatrixXd V_B2;
	Eigen::MatrixXi F_B2;
	igl::readOBJ("ball.obj", V_B2, F_B2);
	const Slicer_2 slicer = current_slicer; //加个const看看有没有涉及更改
	/*slicer.load((file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(cont_number_of_queue) + ".obj").c_str());*/

	/////////////////Get sample points//////////////////
	vector<cv::Point3d> all_sample_points_in_triangles(need_detect_triangle.size());
	for (int i = 0; i < need_detect_triangle.size(); i++) {

		cv::Point3d Point1(slicer.positions[need_detect_triangle[i][0]][0], slicer.positions[need_detect_triangle[i][0]][1], slicer.positions[need_detect_triangle[i][0]][2]);
		cv::Point3d Point2(slicer.positions[need_detect_triangle[i][1]][0], slicer.positions[need_detect_triangle[i][1]][1], slicer.positions[need_detect_triangle[i][1]][2]);
		cv::Point3d Point3(slicer.positions[need_detect_triangle[i][2]][0], slicer.positions[need_detect_triangle[i][2]][1], slicer.positions[need_detect_triangle[i][2]][2]);
		double a = distance3d(Point1, Point2);
		double b = distance3d(Point1, Point3);
		double c = distance3d(Point2, Point3);
		cv::Point3d V_incentre(((a * Point1.x + b * Point2.x + c * Point3.x) / (a + b + c)), ((a * Point1.y + b * Point2.y + c * Point3.y) / (a + b + c)), ((a * Point1.z + b * Point2.z + c * Point3.z) / (a + b + c)));
		all_sample_points_in_triangles[i] = V_incentre;
	}
	///////////////////////////////////////////////////
	cout << "% " << slicer.triangles.size() << endl;
	cout << "% " << need_detect_triangle.size() << endl;
	vector<vector<int>> accessible_ori_of_need_detect_V(need_detect_triangle.size(), vector<int>(sampling_subtractive.sample_points.size(), 10000000));
	vector<vector<Eigen::MatrixXd>> vis_points(1);
	vector<vector<vector<Eigen::Vector3d>>> vis_lines(1);


	cout << "%sampling_subtractive.sample_points.size() " << sampling_subtractive.sample_points.size() << endl;
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {
		//先去除水平面以下的方向  //注意只有最底下的block需要限制!!!!!!!!!!!!!!
		/*if (sampling_subtractive.sample_points[ori].z < 0.2)
			continue;*/
		vector<std::vector<Eigen::Vector3d>> temp_new_V_remain(slicer.triangles.size());
		for (int i = 0; i < slicer.triangles.size(); i++) {
			temp_new_V_remain[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_new_V_remain[i][k](0, 0) = slicer.positions[slicer.triangles[i][k]][0];
				temp_new_V_remain[i][k](1, 0) = slicer.positions[slicer.triangles[i][k]][1];
				temp_new_V_remain[i][k](2, 0) = slicer.positions[slicer.triangles[i][k]][2];
			}
		}

		std::vector<Eigen::Vector3d> temp_V_need_detect(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect[i](0, 0) = all_sample_points_in_triangles[i].x;
			temp_V_need_detect[i](1, 0) = all_sample_points_in_triangles[i].y;
			temp_V_need_detect[i](2, 0) = all_sample_points_in_triangles[i].z;
		}

		vector<std::vector<Eigen::Vector3d>> temp_V_need_detect_triangle(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect_triangle[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_V_need_detect_triangle[i][k](0, 0) = slicer.positions[need_detect_triangle[i][k]][0];
				temp_V_need_detect_triangle[i][k](1, 0) = slicer.positions[need_detect_triangle[i][k]][1];
				temp_V_need_detect_triangle[i][k](2, 0) = slicer.positions[need_detect_triangle[i][k]][2];
			}
		}

		Eigen::Matrix3d rotMatrix;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < temp_new_V_remain.size(); i++)
			for (int j = 0; j < temp_new_V_remain[i].size(); j++)
				temp_new_V_remain[i][j] = rotMatrix.inverse() * temp_new_V_remain[i][j];
		for (int i = 0; i < temp_V_need_detect.size(); i++)
			temp_V_need_detect[i] = rotMatrix.inverse() * temp_V_need_detect[i];
		for (int i = 0; i < temp_V_need_detect_triangle.size(); i++)
			for (int j = 0; j < temp_V_need_detect_triangle[i].size(); j++)
				temp_V_need_detect_triangle[i][j] = rotMatrix.inverse() * temp_V_need_detect_triangle[i][j];

		vector<Eigen::Vector3d> all_normals_of_need_detect_triangle;
		all_normals_of_need_detect_triangle.clear();
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			cv::Point3d Point1(temp_V_need_detect_triangle[i][0](0, 0), temp_V_need_detect_triangle[i][0](1, 0), temp_V_need_detect_triangle[i][0](2, 0));
			cv::Point3d Point2(temp_V_need_detect_triangle[i][1](0, 0), temp_V_need_detect_triangle[i][1](1, 0), temp_V_need_detect_triangle[i][1](2, 0));
			cv::Point3d Point3(temp_V_need_detect_triangle[i][2](0, 0), temp_V_need_detect_triangle[i][2](1, 0), temp_V_need_detect_triangle[i][2](2, 0));
			double na = (Point2.y - Point1.y) * (Point3.z - Point1.z) - (Point2.z - Point1.z) * (Point3.y - Point1.y);
			double nb = (Point2.z - Point1.z) * (Point3.x - Point1.x) - (Point2.x - Point1.x) * (Point3.z - Point1.z);
			double nc = (Point2.x - Point1.x) * (Point3.y - Point1.y) - (Point2.y - Point1.y) * (Point3.x - Point1.x);
			Eigen::Vector3d vn(na, nb, nc);
			vn.normalize();
			all_normals_of_need_detect_triangle.push_back(vn);
			//cout << vn.x() << " "<<vn.y()<<" "<<vn.z() << endl;
		}

		vector<double> max_z_of_triangles(temp_new_V_remain.size());
		for (int i = 0; i < temp_new_V_remain.size(); i++) {
			max_z_of_triangles[i] = MIN_D;
			for (int j = 0; j < 3; j++)
				max_z_of_triangles[i] = max(max_z_of_triangles[i], temp_new_V_remain[i][j](2, 0));
		}

		///////////////////collision detection////////////////////////
		PrepareToolForCollision(cutting_tool);

		for (int i = 0; i < temp_V_need_detect.size(); i++) {
			bool jud_collision = false;
			Eigen::Vector3d center_point;

			center_point.x() = temp_V_need_detect[i](0, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].x();
			center_point.y() = temp_V_need_detect[i](1, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].y();
			center_point.z() = temp_V_need_detect[i](2, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].z();

			for (int ii = 0; ii < temp_new_V_remain.size(); ii++) {

				if (CheckToolCollisionWithCell(center_point, temp_new_V_remain[ii], max_z_of_triangles[ii], cutting_tool, 30.0, 3.0)) {
					jud_collision = true;
					break;
				}
			}

			if (jud_collision == false) {
				accessible_ori_of_need_detect_V[i][ori] = 0;
			}
		}


		//visualize
		/*rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorAfter,vectorBefore).toRotationMatrix();
		for (int i = 0; i < vis_points.size(); i++)
			vis_points[i] = rotMatrix.inverse() * vis_points[i];
		ofstream all_balls(".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + ".obj");
		for (int i = 0; i < vis_points.size(); i++) {
			for (int j = 0; j < V_B2.rows(); j++)
				all_balls << "v " << V_B2(j, 0) + vis_points[i](0, 0) << " " << V_B2(j, 1) + vis_points[i](1, 0) << " " << V_B2(j, 2) + vis_points[i](2, 0) << " 0.9" << " 0.05" << " 0.05" << endl;
			for (int j = 0; j < F_B2.rows(); j++)
				all_balls << "f " << F_B2(j, 0) + i * V_B2.rows() + 1 << " " << F_B2(j, 1) + i * V_B2.rows() + 1 << " " << F_B2(j, 2) + i * V_B2.rows() + 1 << endl;
		}
		all_balls.close();
		Visual vis;
		vis.generateModelForRendering_7(sampling_subtractive.sample_points[ori], ".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + "_orientation.obj");*/
	}


	//////////////////////////////////graph cut////////////////////////////////////////
	//不应该出现cont_ori == 0的情况
	int cont_revise = 0;
	for (int i = 0; i < accessible_ori_of_need_detect_V.size(); i++) {
		int cont_ori = 0;
		for (int j = 0; j < accessible_ori_of_need_detect_V[i].size(); j++) {
			if (accessible_ori_of_need_detect_V[i][j] == 0)
				cont_ori++;
		}
		//暂时先强制任意方向可达
		if (cont_ori == 0) {
			cont_revise++;
			//cout << "*** " << need_detect_triangle[i][0] << endl;
			for (int j = 0; j < accessible_ori_of_need_detect_V[i].size(); j++)
				accessible_ori_of_need_detect_V[i][j] = 0;
		}
		/*if (accessible_ori_of_need_detect_V[i][123] != 0)
			cout << "no" << endl;*/
	}
	cout << cont_revise << endl;


	vector<std::vector<int>> pixels_relations;
	vector<int> length_edges;
	for (int i = 0; i < need_detect_triangle.size(); i++)
		for (int j = 0; j < need_detect_triangle.size(); j++) {
			bool jud_adjacent = false;
			vector<int> temp_edge(2);
			for (int ii = 0; ii < 3; ii++) {
				for (int jj = 0; jj < 3; jj++) {
					if (need_detect_triangle[i][ii] == need_detect_triangle[j][jj]) {
						/*temp_edge[0] = i;
						temp_edge[1] = j;
						pixels_relations.push_back(temp_edge);
						break;*/
						if (need_detect_triangle[i][(ii + 1) % 3] == need_detect_triangle[j][(jj + 1) % 3] || need_detect_triangle[i][(ii + 1) % 3] == need_detect_triangle[j][(jj + 2) % 3]) {
							temp_edge[0] = i;
							temp_edge[1] = j;
							pixels_relations.push_back(temp_edge);
							float temp_length = distanceVec3(slicer.positions[need_detect_triangle[i][ii]], slicer.positions[need_detect_triangle[i][(ii + 1) % 3]]) * 100;
							length_edges.push_back(temp_length);
							jud_adjacent = true;
							break;
						}
						else if (need_detect_triangle[i][(ii + 2) % 3] == need_detect_triangle[j][(jj + 1) % 3] || need_detect_triangle[i][(ii + 2) % 3] == need_detect_triangle[j][(jj + 2) % 3]) {
							temp_edge[0] = i;
							temp_edge[1] = j;
							pixels_relations.push_back(temp_edge);
							float temp_length = distanceVec3(slicer.positions[need_detect_triangle[i][ii]], slicer.positions[need_detect_triangle[i][(ii + 2) % 3]]) * 100;
							length_edges.push_back(temp_length);
							jud_adjacent = true;
							break;
						}
						else {
							jud_adjacent = true;
							break;
						}
					}
				}
				if (jud_adjacent == true)
					break;
			}
			/*if (jud_adjacent == false) {
				temp_edge[0] = i;
				temp_edge[1] = j;
				pixels_relations.push_back(temp_edge);
				length_edges.push_back(100000);
			}*/
		}


	cout << endl << "Graph Cut................" << endl;
	cout << "% " << need_detect_triangle.size() << endl;
	cout << "% " << sampling_subtractive.sample_points.size() << endl;
	vector<int> result = GeneralGraph_DArraySArraySpatVarying(need_detect_triangle.size(), sampling_subtractive.sample_points.size(), accessible_ori_of_need_detect_V, pixels_relations, length_edges);

	for (int i = 0; i < need_detect_triangle.size(); i++)
		for (int j = i + 1; j < need_detect_triangle.size(); j++)
			if (result[i] > result[j]) {
				swap(result[i], result[j]);
				swap(need_detect_triangle[i], need_detect_triangle[j]);
			}

	//cout << "normal:" << result[0]<<endl;
	ofstream ofile(".\\vis\\normal_of_pathches.txt");
	vector<vector<vasco::core::Tri3>> vis_triangles(1);
	vector<vasco::core::Vec3> vis_positions = slicer.positions;

	int cont_patch = 0;
	vector<Eigen::Vector3d> points_in_cell;
	vector<Eigen::Vector3d> normals;
	points_in_cell.push_back(Eigen::Vector3d(slicer.positions[need_detect_triangle[0][0]][0], slicer.positions[need_detect_triangle[0][0]][1], slicer.positions[need_detect_triangle[0][0]][2])); //point_in_cell初始放入第一个点
	normals.push_back(sampling_subtractive.sample_points[result[0]]);
	ofile << normals[0].x() << " " << normals[0].y() << " " << normals[0].z() << endl;


	for (int i = 0; i < need_detect_triangle.size(); i++) {
		if (i != 0 && result[i - 1] < result[i]) {
			cont_patch++;
			vector<Eigen::MatrixXd> temp_vecc(0);
			vis_points.push_back(temp_vecc);
			vector<vector<Eigen::Vector3d>> temp_vec(0);
			vis_lines.push_back(temp_vec);
			vector<vasco::core::Tri3> temp_tri(0);
			vis_triangles.push_back(temp_tri);
			points_in_cell.push_back(Eigen::Vector3d(slicer.positions[need_detect_triangle[i][0]][0], slicer.positions[need_detect_triangle[i][0]][1], slicer.positions[need_detect_triangle[i][0]][2]));
			normals.push_back(sampling_subtractive.sample_points[result[i]]);
			ofile << normals[normals.size() - 1].x() << " " << normals[normals.size() - 1].y() << " " << normals[normals.size() - 1].z() << endl;
		}

		//Eigen::MatrixXd temp_point(3, 1);
		//temp_point(0, 0) = V(index_V_need_detect[i], 0);
		//temp_point(1, 0) = V(index_V_need_detect[i], 1);
		//temp_point(2, 0) = V(index_V_need_detect[i], 2);
		////cout << vis_points[cont_patch].size() << endl;
		//vis_points[cont_patch].push_back(temp_point);
		//cout << "t" << endl;
		vector<Eigen::Vector3d> temp_vec;

		temp_vec.clear();
		vasco::core::Tri3 temp_tri;
		//cout << index_V_need_detect[i] << endl;
		for (int j = 0; j < 3; j++) {
			//cout << "f" << endl;
			Eigen::Vector3d temp_vec_2(slicer.positions[need_detect_triangle[i][j]][0], slicer.positions[need_detect_triangle[i][j]][1], slicer.positions[need_detect_triangle[i][j]][2]);
			//cout << "g" << endl;
			temp_vec.push_back(temp_vec_2);
			temp_tri[j] = need_detect_triangle[i][j];
		}
		//cout << "d" << endl;
		//cout << vis_lines.size() << " " << cont_patch << endl;
		vis_lines[cont_patch].push_back(temp_vec);
		vis_triangles[cont_patch].push_back(temp_tri);
		//cout << vis_lines[cont_patch].size() << endl;
		//cout << "e" << endl;
	}

	//visualize	
	/*ofstream all_balls(".\\vis\\coral_subtractive_decompose.obj");
	int cont_v = 0;
	for (int i = 0; i < vis_points.size(); i++) {
		double r = rand() / double(RAND_MAX);
		double g = rand() / double(RAND_MAX);
		double b = rand() / double(RAND_MAX);
		for (int j = 0; j < vis_points[i].size(); j++) {
			for (int k = 0; k < V_B2.rows(); k++)
				all_balls << "v " << V_B2(k, 0) + vis_points[i][j](0, 0) << " " << V_B2(k, 1) + vis_points[i][j](1, 0) << " " << V_B2(k, 2) + vis_points[i][j](2, 0) << " " << r << " " << g << " " << b << endl;
			for (int k = 0; k < F_B2.rows(); k++)
				all_balls << "f " << F_B2(k, 0) + cont_v * V_B2.rows() + 1 << " " << F_B2(k, 1) + cont_v * V_B2.rows() + 1 << " " << F_B2(k, 2) + cont_v * V_B2.rows() + 1 << endl;
			cont_v++;
		}
	}
	all_balls.close();*/
	for (int t = 0; t < vis_lines.size(); t++) {
		std::ofstream dstream(".\\vis\\patch-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(t) + ".stl");
		if (!dstream.is_open()) {
			std::cout << "can not open " << std::endl;
			return;
		}
		dstream << "solid STL generated by MeshLab" << std::endl;
		for (int i = 0; i < vis_lines[t].size(); i++) {
			dstream << "  facet normal " << "0 0 0" << std::endl;
			dstream << "    outer loop" << std::endl;
			for (int j = 0; j < 3; j++) {
				dstream << "      vertex  " << vis_lines[t][i][j][0] << " " << vis_lines[t][i][j][1] << " " << vis_lines[t][i][j][2] << std::endl;
			}
			dstream << "    endloop" << std::endl;
			dstream << "  endfacet" << std::endl;
		}
		dstream << "endsolid vcg" << std::endl;
		dstream.close();
	}
	for (int t = 0; t < vis_triangles.size(); t++) {
		std::string mesh_name = "subtractive_patch_" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(t);
		polyscope::registerSurfaceMesh(mesh_name, vis_positions, vis_triangles[t]);
	}
	polyscope::show();
	ofile.close();


	cout << cont_patch << endl;
	vector<vector<double>> colors(vis_lines.size(), vector<double>(3));
	Visual vis;
	vis.generateModelForRendering_10(vis_lines, colors, ".\\vis\\coral_subtractive_patch_voronoi_cell.obj");
	Visual vis_2;
	vis_2.generateModelForRendering_11(points_in_cell, normals, colors, ".\\vis\\coral_subtractive_patch_normal.obj");
	//Visual vis;
	//vis.generateModelForRendering_7(sampling_subtractive.sample_points[ori], ".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + "_orientation.obj");
	//////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////
}

void HybridManufacturing::subtractive_accessibility_decomposition_global(int height_of_beam_search, cutter cutting_tool)
{
	sampling_subtractive.OrientationSamplePoints();	//sampling_subtractive生成球面采样点
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 27;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 23;
	cutting_tool.carriage_height = 33;

	/////show all accessible points in every orientation/////
	Eigen::MatrixXd V_B2;
	Eigen::MatrixXi F_B2;
	igl::readOBJ("ball.obj", V_B2, F_B2);

	Slicer_2 slicer_load_current_patch;

	Slicer_2 slicer_load_next_patch;

	// 读取并拼接 .\vis\block_patch-1_.obj ~ .\vis\block_patch-height_of_beam_search_.obj
	std::vector<int> merged_vertex_source_patch_id;
	std::vector<int> merged_face_source_patch_id;

	Slicer_2 slicer_load_merged_patch = MergeBlockPatchesWithDedup(
		height_of_beam_search,
		merged_vertex_source_patch_id,
		merged_face_source_patch_id,
		1e-3);
	slicer_load_merged_patch.save(".\\vis\\block_patch_merged_removedup.obj");
	std::cout << "[Info] merged patch mesh: V=" << slicer_load_merged_patch.positions.size()
		<< ", F=" << slicer_load_merged_patch.triangles.size() << std::endl;

	// 计算 merged patch 上每个三角面在每个采样方向下的最大碰撞 patch_index
	const auto merged_face_min_collision_patch = EvaluateMergedPatchToolCollision(
		slicer_load_merged_patch,
		merged_face_source_patch_id,
		cutting_tool);

	std::cout << "[Info] merged-face min-collision-patch matrix computed: faces="
		<< merged_face_min_collision_patch.size() << ", orientations="
		<< sampling_subtractive.sample_points.size() << std::endl;

	// ---------------- graph cut on merged patch ----------------
	{
		const int face_count = static_cast<int>(slicer_load_merged_patch.triangles.size());
		const int ori_count = static_cast<int>(sampling_subtractive.sample_points.size());
		const int block_count = height_of_beam_search;
		const int num_labels = block_count * ori_count;
		const int INF_COST = 10000000;

		if (face_count == 0 || ori_count == 0 || block_count <= 0) {
			std::cout << "[Warn] skip graph cut: invalid sizes. face_count=" << face_count
				<< ", ori_count=" << ori_count << ", block_count=" << block_count << std::endl;
		}
		else {
			auto encode_label = [ori_count](int patch_id, int ori_id) -> int {
				// patch_id: 1..block_count, ori_id: 0..ori_count-1
				return (patch_id - 1) * ori_count + ori_id;
				};

			auto decode_label = [ori_count](int label_id, int& patch_id, int& ori_id) {
				patch_id = label_id / ori_count + 1;
				ori_id = label_id % ori_count;
				};

			// data cost: face x label
			std::vector<std::vector<int>> data_value(face_count, std::vector<int>(num_labels, INF_COST));

			int warn_a_lt_b = 0;
			int warn_no_feasible = 0;

			for (int i = 0; i < face_count; ++i) {
				int b = 1;
				if (i < static_cast<int>(merged_face_source_patch_id.size())) {
					b = merged_face_source_patch_id[i];
				}
				else {
					std::cout << "[what?]i >= static_cast<int>(merged_face_source_patch_id.size())" << std::endl;
				}
				b = std::max(1, std::min(block_count, b));

				int feasible_cnt = 0;

				for (int ori = 0; ori < ori_count; ++ori) {
					int a = merged_face_min_collision_patch[i][ori];
					if (a == -1) {
						a = 1;
					}

					a = std::max(1, std::min(a, block_count + 1));

					// 按 [a..b] 添加标签，若 a>b 不添加并提示
					if (a > b) {
						++warn_a_lt_b;
						continue;
					}

					for (int p = a; p <= b; ++p) {
						const int lid = encode_label(p, ori);
						data_value[i][lid] = 0;
						++feasible_cnt;
					}
				}

				// 防御：若该面没有任何可行label，给一个保底label避免graph cut退化
				if (feasible_cnt == 0) {
					++warn_no_feasible;
					data_value[i][encode_label(b, 30)] = 0;
				}
			}

			if (warn_a_lt_b > 0) {
				std::cout << "[Warn] graph-cut label range invalid (a<b) count = " << warn_a_lt_b << std::endl;
			}
			if (warn_no_feasible > 0) {
				std::cout << "[Warn] faces with no feasible label = " << warn_no_feasible
					<< " (fallback label assigned)." << std::endl;
			}

			// 建立面邻接图（共享边）
			std::vector<std::vector<int>> pixels_relations;
			std::vector<int> length_edges;

			std::map<std::pair<int, int>, int> edge_owner; // edge -> face id
			edge_owner.clear();

			for (int i = 0; i < face_count; ++i) {
				const auto& tri = slicer_load_merged_patch.triangles[i];
				for (int e = 0; e < 3; ++e) {
					int u = tri[e];
					int v = tri[(e + 1) % 3];
					if (u > v) std::swap(u, v);

					const std::pair<int, int> key(u, v);
					auto it = edge_owner.find(key);
					if (it == edge_owner.end()) {
						edge_owner.insert({ key, i });
					}
					else {
						const int j = it->second;
						if (j != i) {
							pixels_relations.push_back({ j, i });

							const auto& p1 = slicer_load_merged_patch.positions[u];
							const auto& p2 = slicer_load_merged_patch.positions[v];
							const int w = std::max(1, static_cast<int>(distanceVec3(p1, p2) * 100.0));
							length_edges.push_back(w);
						}
					}
				}
			}

			std::cout << "[Info] graph-cut nodes=" << face_count
				<< ", labels=" << num_labels
				<< ", edges=" << pixels_relations.size() << std::endl;

			// graph cut
			std::vector<int> result_labels = GeneralGraph_DArraySArraySpatVarying(
				face_count,
				num_labels,
				data_value,
				pixels_relations,
				length_edges);

			// 可视化：按最终label所属patch着色（离散固定色）
			const std::string gc_obj_file = ".\\vis\\merged_patch_graphcut_label.obj";
			std::ofstream ofs(gc_obj_file);
			if (!ofs.is_open()) {
				std::cout << "[Warn] cannot open file for writing: " << gc_obj_file << std::endl;
			}
			else {
				auto color_from_patch = [](int ori_id) -> std::array<double, 3> {
					static const std::array<std::array<double, 3>, 32> palette = { {
						{0.894, 0.102, 0.110}, {0.216, 0.494, 0.722}, {0.302, 0.686, 0.290}, {0.596, 0.306, 0.639},
						{1.000, 0.498, 0.000}, {1.000, 1.000, 0.200}, {0.651, 0.337, 0.157}, {0.969, 0.506, 0.749},
						{0.600, 0.600, 0.600}, {0.121, 0.466, 0.705}, {0.682, 0.780, 0.909}, {1.000, 0.733, 0.470},
						{0.172, 0.627, 0.172}, {0.839, 0.153, 0.157}, {0.580, 0.404, 0.741}, {0.549, 0.337, 0.294},
						{0.890, 0.467, 0.761}, {0.498, 0.498, 0.498}, {0.737, 0.741, 0.133}, {0.090, 0.745, 0.811},
						{0.400, 0.760, 0.647}, {0.988, 0.553, 0.384}, {0.553, 0.627, 0.796}, {0.906, 0.541, 0.765},
						{0.651, 0.847, 0.329}, {1.000, 0.851, 0.184}, {0.898, 0.768, 0.580}, {0.702, 0.702, 0.702},
						{0.984, 0.603, 0.600}, {0.800, 0.922, 0.773}, {0.871, 0.796, 0.894}, {0.996, 0.851, 0.651}
					} };

					int idx = ori_id;
					if (idx < static_cast<int>(palette.size())) {
						return palette[static_cast<size_t>(idx)];
					}

					return palette[idx % palette.size()];

					};

				// 同时输出 face -> (patch,ori)
				std::ofstream lfs(".\\vis\\merged_patch_graphcut_label.txt");
				lfs << "face_id patch_id ori_id label_id\n";

				int v_count = 0;
				int skipped_faces = 0;
				int invalid_vertices = 0;

				// 按 patch_id 收集面，供 polyscope 分 mesh 可视化
				std::map<int, std::vector<vasco::core::Tri3>> patch_to_faces;

				for (int i = 0; i < face_count; ++i) {
					int patch_id = 1, ori_id = 0;
					decode_label(result_labels[i], patch_id, ori_id);

					// 修正：按 patch_id 取颜色
					const auto color = color_from_patch(result_labels[i]);
					const auto& tri = slicer_load_merged_patch.triangles[i];

					// tri 索引与坐标合法性检查
					bool face_ok = true;
					for (int k = 0; k < 3; ++k) {
						const int vid = tri[k];
						if (vid < 0 || vid >= static_cast<int>(slicer_load_merged_patch.positions.size())) {
							face_ok = false;
							break;
						}
						const auto& p = slicer_load_merged_patch.positions[vid];
						if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
							face_ok = false;
							++invalid_vertices;
							break;
						}
					}

					if (!face_ok) {
						++skipped_faces;
						continue;
					}

					for (int k = 0; k < 3; ++k) {
						const auto& p = slicer_load_merged_patch.positions[tri[k]];
						ofs << "v " << p[0] << " " << p[1] << " " << p[2] << " "
							<< color[0] << " " << color[1] << " " << color[2] << "\n";
					}
					ofs << "f " << (v_count + 1) << " " << (v_count + 2) << " " << (v_count + 3) << "\n";
					v_count += 3;

					lfs << i << " " << patch_id << " " << ori_id << " " << result_labels[i] << "\n";

					// 保存到 patch 子网格
					patch_to_faces[result_labels[i]].push_back(tri);
				}

				if (skipped_faces > 0) {
					std::cout << "[Warn] skipped faces while writing OBJ: " << skipped_faces
						<< ", invalid_vertices=" << invalid_vertices << std::endl;
				}

				// polyscope: 按 patch_id 拆分 mesh，并设置固定颜色
				std::vector<Eigen::Vector3d> patch_arrow_points;
				std::vector<Eigen::Vector3d> patch_arrow_dirs;
				std::vector<std::vector<double>> patch_arrow_colors;

				Eigen::Vector3d bb_min(MAX_D, MAX_D, MAX_D), bb_max(MIN_D, MIN_D, MIN_D);
				for (const auto& p : slicer_load_merged_patch.positions) {
					bb_min.x() = std::min(bb_min.x(), p[0]);
					bb_min.y() = std::min(bb_min.y(), p[1]);
					bb_min.z() = std::min(bb_min.z(), p[2]);
					bb_max.x() = std::max(bb_max.x(), p[0]);
					bb_max.y() = std::max(bb_max.y(), p[1]);
					bb_max.z() = std::max(bb_max.z(), p[2]);
				}
				const double model_diag = (bb_max - bb_min).norm();
				const double arrow_start_offset = std::max(1.0, 0.06 * model_diag);
				const double arrow_length = std::max(2.0, 0.12 * model_diag);

				for (const auto& kv : patch_to_faces) {
					const int patch_id = kv.first;
					const auto& patch_faces = kv.second;
					if (patch_faces.empty()) continue;
					int de_patch;
					int ori_id = 0;
					decode_label(patch_id, de_patch, ori_id);
					const std::string mesh_name =
						"merged_patch_graphcut_patch_" + std::to_string(height_of_beam_search) + "_" + std::to_string(de_patch) + "_" + std::to_string(ori_id);

					auto* ps_mesh = polyscope::registerSurfaceMesh(
						mesh_name,
						slicer_load_merged_patch.positions,
						patch_faces);

					const auto c = color_from_patch(ori_id);
					ps_mesh->setSurfaceColor({ (float)c[0], (float)c[1], (float)c[2] });

					std::cout << "[debug]mesh_name = " << mesh_name << ", patch_id=" << patch_id << ", de_patch=" << de_patch << ", ori_id=" << ori_id
						<< ", c= " << c[0] << "," << c[1] << "," << c[2] << std::endl;

					Eigen::Vector3d patch_center(0.0, 0.0, 0.0);
					int center_count = 0;
					for (const auto& tri : patch_faces) {
						for (int k = 0; k < 3; ++k) {
							const int vid = tri[k];
							if (vid >= 0 && vid < static_cast<int>(slicer_load_merged_patch.positions.size())) {
								const auto& p = slicer_load_merged_patch.positions[vid];
								patch_center.x() += p[0];
								patch_center.y() += p[1];
								patch_center.z() += p[2];
								++center_count;
							}
						}
					}
					if (center_count > 0) {
						patch_center /= static_cast<double>(center_count);

						Eigen::Vector3d ori_dir(0.0, 0.0, 1.0);
						if (ori_id >= 0 && ori_id < static_cast<int>(sampling_subtractive.sample_points.size())) {
							ori_dir = sampling_subtractive.sample_points[ori_id];
						}
						Eigen::Vector3d ori_vec(ori_dir);
						double n = ori_vec.norm();
						if (n < 1e-12) {
							ori_vec = Eigen::Vector3d(0.0, 0.0, 1.0);
						}
						else {
							ori_vec /= n;
						}

						const Eigen::Vector3d arrow_start = patch_center + ori_vec * arrow_start_offset;
						const Eigen::Vector3d arrow_end = arrow_start + ori_vec * arrow_length;

						patch_arrow_points.push_back(arrow_start);
						patch_arrow_dirs.push_back(ori_vec);

						patch_arrow_colors.push_back({ c[0], c[1], c[2] });

						const std::string arrow_name =
							"merged_patch_graphcut_arrow_" + std::to_string(height_of_beam_search) + "_" + std::to_string(de_patch) + "_" + std::to_string(ori_id);

						Eigen::Vector3d ref_axis = (std::abs(ori_vec.z()) < 0.9)
							? Eigen::Vector3d(0.0, 0.0, 1.0)
							: Eigen::Vector3d(0.0, 1.0, 0.0);
						Eigen::Vector3d perp1 = ori_vec.cross(ref_axis);
						double perp1_norm = perp1.norm();
						if (perp1_norm < 1e-12) {
							perp1 = Eigen::Vector3d(1.0, 0.0, 0.0);
						}
						else {
							perp1 /= perp1_norm;
						}

						const double shaft_w = std::max(0.2, arrow_length * 0.08);
						const double head_len = arrow_length * 0.35;
						const double head_w = shaft_w * 2.0;

						Eigen::Vector3d m = arrow_end - ori_vec * head_len;
						Eigen::Vector3d s0 = arrow_start + perp1 * shaft_w;
						Eigen::Vector3d s1 = arrow_start - perp1 * shaft_w;
						Eigen::Vector3d m0 = m + perp1 * shaft_w;
						Eigen::Vector3d m1 = m - perp1 * shaft_w;
						Eigen::Vector3d h0 = m + perp1 * head_w;
						Eigen::Vector3d h1 = m - perp1 * head_w;

						std::vector<vasco::core::Vec3> arrow_pos = {
							{ s0.x(), s0.y(), s0.z() },
							{ s1.x(), s1.y(), s1.z() },
							{ m1.x(), m1.y(), m1.z() },
							{ m0.x(), m0.y(), m0.z() },
							{ h0.x(), h0.y(), h0.z() },
							{ h1.x(), h1.y(), h1.z() },
							{ arrow_end.x(), arrow_end.y(), arrow_end.z() }
						};
						std::vector<vasco::core::Tri3> arrow_tri = {
							{ 0, 1, 2 }, { 0, 2, 3 },
							{ 4, 5, 6 }, { 3, 4, 6 }, { 2, 6, 5 }
						};

						auto* ps_arrow = polyscope::registerSurfaceMesh(arrow_name, arrow_pos, arrow_tri);
						ps_arrow->setSurfaceColor({ (float)c[0], (float)c[1], (float)c[2] });
					}
				}

				if (!patch_arrow_points.empty()) {
					Visual vis_patch_ori;
					vis_patch_ori.generateModelForRendering_11(
						patch_arrow_points,
						patch_arrow_dirs,
						patch_arrow_colors,
						".\\vis\\merged_patch_graphcut_patch_ori_arrows.obj");
				}
				polyscope::show();

				std::cout << "[Info] wrote graph-cut colored OBJ: " << gc_obj_file << std::endl;
				std::cout << "[Info] wrote graph-cut labels txt: .\\vis\\merged_patch_graphcut_label.txt" << std::endl;
			}
		}
	}

}

int HybridManufacturing::subtractive_accessibility_decomposition_local(int height_of_beam_search, cutter cutting_tool)
{
	if (sampling_subtractive.sample_points.empty()) {
		sampling_subtractive.OrientationSamplePoints();	//sampling_subtractive生成球面采样点
	}
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 27;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 23;
	cutting_tool.carriage_height = 33;

	/////show all accessible points in every orientation/////
	Eigen::MatrixXd V_B2;
	Eigen::MatrixXi F_B2;
	igl::readOBJ("ball.obj", V_B2, F_B2);

	Slicer_2 slicer_load_current_patch;

	Slicer_2 slicer_load_next_patch;

	// 读取并拼接 .\vis\block_patch-1_.obj ~ .\vis\block_patch-height_of_beam_search_.obj
	std::vector<int> merged_vertex_source_patch_id;
	std::vector<int> merged_face_source_patch_id;

	Slicer_2 slicer_load_merged_patch = MergeBlockPatchesWithDedup(
		height_of_beam_search,
		merged_vertex_source_patch_id,
		merged_face_source_patch_id,
		1e-5);
	slicer_load_merged_patch.save(".\\vis\\block_patch_merged_removedup_local.obj");
	std::cout << "[Info] merged patch mesh: V=" << slicer_load_merged_patch.positions.size()
		<< ", F=" << slicer_load_merged_patch.triangles.size() << std::endl;

	// 计算 merged patch 上每个三角面在每个采样方向下的最大碰撞 patch_index
	const auto merged_face_min_collision_patch = EvaluateMergedPatchToolCollision(
		slicer_load_merged_patch,
		merged_face_source_patch_id,
		cutting_tool,
		true);

	std::cout << "[Info] merged-face min-collision-patch matrix computed: faces="
		<< merged_face_min_collision_patch.size() << ", orientations="
		<< sampling_subtractive.sample_points.size() << std::endl;

	int ret = 999;
	// ---------------- graph cut on merged patch ----------------
	{
		const int face_count = static_cast<int>(slicer_load_merged_patch.triangles.size());
		const int ori_count = static_cast<int>(sampling_subtractive.sample_points.size());
		const int block_count = height_of_beam_search;
		const int num_labels = block_count * ori_count;
		const int INF_COST = 10000000;

		if (face_count == 0 || ori_count == 0 || block_count <= 0) {
			std::cout << "[Warn] skip graph cut: invalid sizes. face_count=" << face_count
				<< ", ori_count=" << ori_count << ", block_count=" << block_count << std::endl;
		}
		else {
			auto encode_label = [ori_count](int patch_id, int ori_id) -> int {
				// patch_id: 1..block_count, ori_id: 0..ori_count-1
				return (patch_id - 1) * ori_count + ori_id;
				};

			auto decode_label = [ori_count](int label_id, int& patch_id, int& ori_id) {
				patch_id = label_id / ori_count + 1;
				ori_id = label_id % ori_count;
				};

			// data cost: face x label
			std::vector<std::vector<int>> data_value(face_count, std::vector<int>(num_labels, INF_COST));

			int warn_a_lt_b = 0;
			int warn_no_feasible = 0;

			for (int i = 0; i < face_count; ++i) {
				int b = 1;
				if (i < static_cast<int>(merged_face_source_patch_id.size())) {
					b = merged_face_source_patch_id[i];
				}
				else {
					std::cout << "[what?]i >= static_cast<int>(merged_face_source_patch_id.size())" << std::endl;
				}
				b = std::max(1, std::min(block_count, b));

				int feasible_cnt = 0;

				for (int ori = 0; ori < ori_count; ++ori) {
					int a = merged_face_min_collision_patch[i][ori];
					if (a == -1) {
						a = 1;
					}

					a = std::max(1, std::min(a, block_count));

					// 按 [a..b] 添加标签，若 a>b 不添加并提示
					if (a > b) {
						++warn_a_lt_b;
						continue;
					}

					for (int p = a; p <= b; ++p) {
						const int lid = encode_label(p, ori);
						data_value[i][lid] = 0;
						++feasible_cnt;
					}
				}

				// 防御：若该面没有任何可行label，给一个保底label避免graph cut退化
				if (feasible_cnt == 0) {
					++warn_no_feasible;
					data_value[i][encode_label(b, 30)] = 0;
				}
			}

			if (warn_a_lt_b > 0) {
				std::cout << "[Warn] graph-cut label range invalid (a<b) count = " << warn_a_lt_b << std::endl;
			}
			if (warn_no_feasible > 0) {
				std::cout << "[Warn] faces with no feasible label = " << warn_no_feasible
					<< " (fallback label assigned)." << std::endl;
			}

			// 建立面邻接图（共享边）
			std::vector<std::vector<int>> pixels_relations;
			std::vector<int> length_edges;

			std::map<std::pair<int, int>, int> edge_owner; // edge -> face id
			edge_owner.clear();

			for (int i = 0; i < face_count; ++i) {
				const auto& tri = slicer_load_merged_patch.triangles[i];
				for (int e = 0; e < 3; ++e) {
					int u = tri[e];
					int v = tri[(e + 1) % 3];
					if (u > v) std::swap(u, v);

					const std::pair<int, int> key(u, v);
					auto it = edge_owner.find(key);
					if (it == edge_owner.end()) {
						edge_owner.insert({ key, i });
					}
					else {
						const int j = it->second;
						if (j != i) {
							pixels_relations.push_back({ j, i });

							const auto& p1 = slicer_load_merged_patch.positions[u];
							const auto& p2 = slicer_load_merged_patch.positions[v];
							const int w = std::max(1, static_cast<int>(distanceVec3(p1, p2) * 100.0));
							length_edges.push_back(w);
						}
					}
				}
			}

			std::cout << "[Info] graph-cut nodes=" << face_count
				<< ", labels=" << num_labels
				<< ", edges=" << pixels_relations.size() << std::endl;

			// graph cut
			std::vector<int> result_labels = GeneralGraph_DArraySArraySpatVarying(
				face_count,
				num_labels,
				data_value,
				pixels_relations,
				length_edges);

			auto label_type_count = result_labels;
			std::sort(label_type_count.begin(), label_type_count.end());
			ret = std::unique(label_type_count.begin(), label_type_count.end()) - label_type_count.begin();

		}
	}
	return ret;
}

vector<vector<int>> HybridManufacturing::getAccessOri(const Slicer_2& slicer, Slicer_2& slicer_load_patch, vector<vasco::core::Vec3>& all_sample_points_in_triangles, cutter cutting_tool)
{
	auto& need_detect_triangle = slicer_load_patch.triangles;
	vector<vector<int>> accessible_ori_of_need_detect_V(need_detect_triangle.size(), vector<int>(sampling_subtractive.sample_points.size(), 10000000));
	cout << "%sampling_subtractive.sample_points.size() " << sampling_subtractive.sample_points.size() << endl;
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {
		vector<std::vector<Eigen::Vector3d>> temp_new_V_remain(slicer.triangles.size());
		for (int i = 0; i < slicer.triangles.size(); i++) {
			temp_new_V_remain[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_new_V_remain[i][k](0, 0) = slicer.positions[slicer.triangles[i][k]][0];
				temp_new_V_remain[i][k](1, 0) = slicer.positions[slicer.triangles[i][k]][1];
				temp_new_V_remain[i][k](2, 0) = slicer.positions[slicer.triangles[i][k]][2];
			}
		}

		std::vector<Eigen::Vector3d> temp_V_need_detect(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect[i](0, 0) = all_sample_points_in_triangles[i][0];
			temp_V_need_detect[i](1, 0) = all_sample_points_in_triangles[i][1];
			temp_V_need_detect[i](2, 0) = all_sample_points_in_triangles[i][2];
		}

		vector<std::vector<Eigen::Vector3d>> temp_V_need_detect_triangle(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect_triangle[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_V_need_detect_triangle[i][k](0, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][0];
				temp_V_need_detect_triangle[i][k](1, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][1];
				temp_V_need_detect_triangle[i][k](2, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][2];
			}
		}

		Eigen::Matrix3d rotMatrix;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < temp_new_V_remain.size(); i++)
			for (int j = 0; j < temp_new_V_remain[i].size(); j++)
				temp_new_V_remain[i][j] = rotMatrix.inverse() * temp_new_V_remain[i][j];
		for (int i = 0; i < temp_V_need_detect.size(); i++)
			temp_V_need_detect[i] = rotMatrix.inverse() * temp_V_need_detect[i];
		for (int i = 0; i < temp_V_need_detect_triangle.size(); i++)
			for (int j = 0; j < temp_V_need_detect_triangle[i].size(); j++)
				temp_V_need_detect_triangle[i][j] = rotMatrix.inverse() * temp_V_need_detect_triangle[i][j];

		vector<Eigen::Vector3d> all_normals_of_need_detect_triangle;
		all_normals_of_need_detect_triangle.clear();
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			cv::Point3d Point1(temp_V_need_detect_triangle[i][0](0, 0), temp_V_need_detect_triangle[i][0](1, 0), temp_V_need_detect_triangle[i][0](2, 0));
			cv::Point3d Point2(temp_V_need_detect_triangle[i][1](0, 0), temp_V_need_detect_triangle[i][1](1, 0), temp_V_need_detect_triangle[i][1](2, 0));
			cv::Point3d Point3(temp_V_need_detect_triangle[i][2](0, 0), temp_V_need_detect_triangle[i][2](1, 0), temp_V_need_detect_triangle[i][2](2, 0));
			double na = (Point2.y - Point1.y) * (Point3.z - Point1.z) - (Point2.z - Point1.z) * (Point3.y - Point1.y);
			double nb = (Point2.z - Point1.z) * (Point3.x - Point1.x) - (Point2.x - Point1.x) * (Point3.z - Point1.z);
			double nc = (Point2.x - Point1.x) * (Point3.y - Point1.y) - (Point2.y - Point1.y) * (Point3.x - Point1.x);
			Eigen::Vector3d vn(na, nb, nc);
			vn.normalize();
			all_normals_of_need_detect_triangle.push_back(vn);
			//cout << vn.x() << " "<<vn.y()<<" "<<vn.z() << endl;
		}

		vector<double> max_z_of_triangles(temp_new_V_remain.size());
		for (int i = 0; i < temp_new_V_remain.size(); i++) {
			max_z_of_triangles[i] = MIN_D;
			for (int j = 0; j < 3; j++)
				max_z_of_triangles[i] = max(max_z_of_triangles[i], temp_new_V_remain[i][j](2, 0));
		}

		///////////////////collision detection////////////////////////
		PrepareToolForCollision(cutting_tool);

		for (int i = 0; i < temp_V_need_detect.size(); i++) {
			Eigen::Vector3d center_point;

			center_point.x() = temp_V_need_detect[i](0, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].x();
			center_point.y() = temp_V_need_detect[i](1, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].y();
			center_point.z() = temp_V_need_detect[i](2, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].z();

			bool jud_collision = false;

			for (int ii = 0; ii < temp_new_V_remain.size(); ii++) {
				if (CheckToolCollisionWithCell(center_point, temp_new_V_remain[ii], max_z_of_triangles[ii], cutting_tool, 30.0, 3.0)) {
					jud_collision = true;
					break;
				}
			}

			if (jud_collision == false) {
				//vis_points.push_back(temp_V_need_detect[i]);
				accessible_ori_of_need_detect_V[i][ori] = 0;
			}
			//else
				//cout << "% " << endl;
		}
	}
	return accessible_ori_of_need_detect_V;
}

void HybridManufacturing::EvaluateCandidateNode(
	int& i,
	vector<int>& candidate_nodes,
	OrientationScores& all_calculated_value,
	vector<BeamTreeEntry>& tree_entries,
	const vector<vector<Eigen::Vector3d>>& save_ori,
	int height_of_beam_search,
	int cont_number_of_queue,
	const string& file_name,
	int now_last_node,
	nozzle the_nozzle)
{
	printf("\r[%d%%]>", i * 100 / (candidate_nodes.size() - 1));
	for (int j = 1; j <= i * 20 / candidate_nodes.size(); j++)
		cout << "▉";

	int index_of_pre_node = tree_entries[candidate_nodes[i]].pre_queue_index;
	vector<vector<cv::Point3d>> all_cut_layers;
	vector<int> all_cut_layers_dependency_layer;
	all_cut_layers.clear();
	all_cut_layers_dependency_layer.clear();
	for (int j = 0; j < tree_entries[candidate_nodes[i]].cut_layers.size(); j++) {
		int index_of_layers = tree_entries[candidate_nodes[i]].cut_layers[j];
		all_cut_layers.push_back(tree_entries[candidate_nodes[i]].layers[index_of_layers]);
		all_cut_layers_dependency_layer.push_back(tree_entries[candidate_nodes[i]].cut_layers_dependency_layer[j]);
	}

	Slicer_2 slicer_G;
	if (tree_entries[now_last_node].judge_continue == false)
		slicer_G.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + ".obj");
	else
		slicer_G.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + "_" + to_string(tree_entries[candidate_nodes[i]].continue_id - 1) + "_subblock.obj");

	all_calculated_value[i] = GainMesh(
		slicer_G,
		all_cut_layers,
		save_ori[candidate_nodes[i]][0],
		height_of_beam_search,
		cont_number_of_queue,
		index_of_pre_node,
		all_cut_layers_dependency_layer,
		tree_entries[now_last_node].judge_continue,
		tree_entries[candidate_nodes[i]].continue_id);

	calculate_fragile_value(all_calculated_value[i], all_cut_layers, tree_entries[candidate_nodes[i]].fragile_v);

	tree_entries[candidate_nodes[i]].larger_base = all_calculated_value[i].large_base; //有点疑问，看看candidate_nodes[i]是否一定等于i
	if (all_calculated_value[i].value_of_self_support == 0) {
		candidate_nodes.erase(candidate_nodes.begin() + i);
		all_calculated_value.erase(all_calculated_value.begin() + i);
		i--;
		return;
	}

	detect_collision_with_printing_platform(
		i,
		candidate_nodes,
		all_calculated_value,
		all_cut_layers,
		save_ori[candidate_nodes[i]][0],
		the_nozzle);
}

void HybridManufacturing::LogSelectedCandidateMetrics(
	int cont_w,
	const vector<int>& candidate_nodes,
	int now_last_node,
	const vector<vector<Eigen::Vector3d>>& save_ori,
	const vector<BeamTreeEntry>& tree_entries,
	queue<int>& last_step_nodes,
	const OrientationScores& all_calculated_value,
	const OrientationScores& pure_value,
	int& sum_candidate_blocks,
	int& sum_connected_components,
	int& cont_extra_additive_orientation,
	vector<double>& evaluation_value)
{
	cout << endl << "selected orientation: " << save_ori[candidate_nodes[cont_w]][0].x() << " " << save_ori[candidate_nodes[cont_w]][0].y() << " " << save_ori[candidate_nodes[cont_w]][0].z();
	cout << endl << "number of candidate_nodes: " << candidate_nodes.size() << endl;
	sum_candidate_blocks += candidate_nodes.size();
	sum_connected_components += tree_entries[candidate_nodes[cont_w]].cut_layers.size();
	if (tree_entries[now_last_node].judge_continue == false)
		last_step_nodes.push(candidate_nodes[cont_w]);
	cout << "self support value of selected node: " << all_calculated_value[cont_w].value_of_self_support << endl << "■■■■■■■■■■■■■■■■■■■" << endl;
	cout << "value_of_more_slice_layers: " << all_calculated_value[cont_w].value_of_more_slice_layers << endl;
	cout << "value_of_area_S: " << all_calculated_value[cont_w].value_of_area_S << endl;
	cout << "value_of_covering_points: " << all_calculated_value[cont_w].value_of_covering_points << endl;
	cout << "value_of_less_clipping_plane: " << all_calculated_value[cont_w].value_of_less_clipping_plane << endl;
	cout << "value_of_fragile: " << all_calculated_value[cont_w].value_of_fragile << endl;
	cout << "value_of_orientation: " << all_calculated_value[cont_w].value_of_orientation << endl;
	cout << "value_of_projected: " << all_calculated_value[cont_w].value_of_projected << endl;
	cout << "■■■■■■■■■■■■■■■■■■■" << endl;

	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	cout << "Pure value orien:" << pure_value[cont_w].value_of_orientation << endl;
	cout << "Pure value fragile:" << pure_value[cont_w].value_of_fragile << endl;
	cout << "Pure value projection:" << pure_value[cont_w].value_of_projected << endl;
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	if (tree_entries[candidate_nodes[cont_w]].judge_continue == true) {
		cout << "*****NEED CHANGE ORIENTATION*****" << endl;
		cont_extra_additive_orientation++;
	}

	evaluation_value[0] += all_calculated_value[cont_w].value_of_more_slice_layers;
	evaluation_value[1] += all_calculated_value[cont_w].value_of_covering_points;
	evaluation_value[2] += all_calculated_value[cont_w].value_of_less_clipping_plane;
	evaluation_value[3] += all_calculated_value[cont_w].value_of_orientation;
	evaluation_value[4] += all_calculated_value[cont_w].value_of_fragile;
	evaluation_value[5] += all_calculated_value[cont_w].value_of_projected;
}

void HybridManufacturing::BuildCutLayersForCandidate(
	int cont_w,
	const vector<int>& candidate_nodes,
	const vector<BeamTreeEntry>& tree_entries,
	vector<vector<cv::Point3d>>& all_cut_layers,
	vector<int>& all_cut_layers_dependency_layer,
	vector<int>& flag_cut_layers_is_hole) const
{
	all_cut_layers.clear();
	all_cut_layers_dependency_layer.clear();
	flag_cut_layers_is_hole.clear();
	for (int j = 0; j < tree_entries[candidate_nodes[cont_w]].cut_layers.size(); j++) {
		int index_of_layers = tree_entries[candidate_nodes[cont_w]].cut_layers[j];
		all_cut_layers.push_back(tree_entries[candidate_nodes[cont_w]].layers[index_of_layers]);
		flag_cut_layers_is_hole.push_back(-1);
		if (tree_entries[candidate_nodes[cont_w]].contain[index_of_layers].size() != 0) {
			all_cut_layers.push_back(tree_entries[candidate_nodes[cont_w]].contain[index_of_layers]);
			all_cut_layers_dependency_layer.push_back(tree_entries[candidate_nodes[cont_w]].cut_layers_dependency_layer[j]);
			flag_cut_layers_is_hole.push_back(all_cut_layers.size() - 2);
		}
		all_cut_layers_dependency_layer.push_back(tree_entries[candidate_nodes[cont_w]].cut_layers_dependency_layer[j]);
	}
}

void HybridManufacturing::PrepareOuterBeamSearchNode(
	queue<int>& last_step_nodes,
	vector<vector<bool>>& is_fragile_V_2,
	int& now_last_node,
	int height_of_beam_search,
	int cont_number_of_queue,
	const vector<BeamTreeEntry>& tree_entries,
	double& sum_time_5)
{
	clock_t start_time_2 = clock();
	is_fragile_V_2 = is_fragile_V;	//fragile信息的快照
	now_last_node = last_step_nodes.front();
	/*if (Tree_nodes_error[now_last_node] == true) {
		last_step_nodes.pop();
		continue;
	}*/

	cout << endl << "Slicing and inner-beam-search......" << endl;
	//load blocks//
	const char* config_path = "config.ini";
	Katana::Instance().config.loadConfig(config_path);
	current_node_mesh.clear();
	if (tree_entries[now_last_node].judge_continue == false) {
		std::string current_file_name = file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl";
		Katana::Instance().stl.loadStl(current_file_name.c_str());	//加载当前节点对应的模型
		LoadAndCheckCurrentNodeMesh(current_file_name, current_node_mesh);	//加载当前节点对应的模型

		cout << "T" << endl;
	}
	else {
		std::string current_file_name = file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(tree_entries[now_last_node].pre_queue_index) + ".stl";
		Katana::Instance().stl.loadStl(current_file_name.c_str());
		LoadAndCheckCurrentNodeMesh(current_file_name, current_node_mesh);
		cout << "U" << now_last_node << endl;
	}
	Katana::Instance().temp_vertices.resize(Katana::Instance().vertices.size());	//更新当前模型的顶点信息
	Katana::Instance().temp_vertices = Katana::Instance().vertices;
	clock_t end_time_2 = clock();
	sum_time_5 += double(end_time_2 - start_time_2) / CLOCKS_PER_SEC;
}

void HybridManufacturing::outer_beam_search(nozzle the_nozzle, cutter cutting_tool)
{
	int W1 = 1;  //4
	clock_t start_time_total, end_time_total;
	clock_t start_time, end_time;
	clock_t start_time_6, end_time_6;
	clock_t start_time_7, end_time_7;
	clock_t start_time_8, end_time_8;
	float total_time_1 = 0, total_time_2 = 0, total_time_3 = 0, total_time_4 = 0;
	double sum_time = 0;
	double sum_time_2 = 0;
	double sum_time_3 = 0;
	double sum_time_4 = 0;
	double sum_time_5 = 0;
	double sum_time_6 = 0;
	int Sum_candidate_blocks = 0;
	int Sum_connected_components = 0;
	vector<double> evaluation_value(6, 0);
	cont_extra_additive_orientation = 0;

	start_time_total = clock();
	bool jud_outer_beam_search_terminate = false;
	vector<BeamTreeEntry> tree_entries;
	vector<int> candidate_nodes;	//当前层生成的候选节点
	queue<int> last_step_nodes;
	vector<vector<Eigen::Vector3d>> save_ori;
	const vector<int> root_cut_layers;
	const vector<Eigen::Vector3d> root_orientations{ Eigen::Vector3d::Zero() };
	BeamTreeEntry root_entry{};
	root_entry.parent_id = -1;
	root_entry.pre_queue_index = -1;
	tree_entries.push_back(root_entry);
	last_step_nodes.push(0);	//推入根节点，作为第一层搜索的初始节点
	tree_entries[0].judge_continue = false;
	tree_entries[0].continue_id = -1;
	tree_entries[0].error = false;
	save_ori.push_back(root_orientations);

	vector<vector<area_S>> ori_all_the_area_S = all_the_area_S;	//不可达点的all_the_area_S的拷贝
	ori_all_the_covering_points = all_the_covering_points;	//不可达点的all_the_covering_points的拷贝
	SAMPLE_ON_BALL sampling;
	int index_node = 0;
	sampling.OrientationSamplePoints_2();	//生成球面采样方向，作为增材打印方向候选
	//Visual Vis_ori;
	//Vis_ori.generateModelForRendering_6(sampling.sample_points, file_name);


	GetALLFragileVertex(sampling);

	Eigen::Matrix3d rotMatrix;
	vector<vector<int>> final_pathes_include_S;
	vector<vector<int>> final_pathes_include_sample_points;
	vector<vector<int>> final_pathes_include_covering_points;
	final_pathes_include_S.push_back(root_cut_layers);
	final_pathes_include_sample_points.push_back(root_cut_layers);
	final_pathes_include_covering_points.push_back(root_cut_layers);
	all_saved_mesh.resize(W1);

	int cont_number_of_queue = 0;	//记录当前处理到第几个节点，用于设定读取的文件名等
	int height_of_beam_search = 0;

	ofstream ofile_cont("E:\\Hybrid manufacturing\\HybridManufacturing\\HybridManufacturing\\models\\coral\\cont_red_number.txt");
	ofstream ofile_cont_size("E:\\Hybrid manufacturing\\HybridManufacturing\\HybridManufacturing\\models\\coral\\cont_size.txt");
	/*bool* flag_sample_point_used = new bool[sampling.sample_points.size()];
	for (int i = 0; i < sampling.sample_points.size(); i++)
		flag_sample_point_used[i] = false;*/
	cout << "--------------------- height of outer-beam search: 1 ---------------------";
	while (last_step_nodes.size() != 0 || candidate_nodes.size() != 0) {
		/*if (height_of_beam_search ==2)
			W1 = 5;
		else
			W1 = 1;*/
		vector<vector<bool>> is_fragile_V_2;
		int now_last_node = -1;
		PrepareOuterBeamSearchNode(
			last_step_nodes,
			is_fragile_V_2,
			now_last_node,
			height_of_beam_search,
			cont_number_of_queue,
			tree_entries,
			sum_time_5);
		//////////////

		vector<bool> judge_S_be_searched;
		vector<bool> judge_covering_points_be_searched;
		sum_time_3 += UpdateSubtractiveDependencyGraph(
			ori_all_the_area_S,
			is_fragile_V_2,
			now_last_node,
			height_of_beam_search,
			cont_number_of_queue,
			tree_entries,
			judge_S_be_searched,
			judge_covering_points_be_searched);
		/////////////////////////////////////

		/*int end;
		if (height_of_beam_search != 4)
			end = sampling.sample_points.size();
		else
			end = 1;*/
		std::cout << sampling.sample_points.size() << endl;
		for (int ori = 0; ori < sampling.sample_points.size(); ori++) {	//枚举所有采样方向，作为增材分层方向
			printf("\r[%d%%]>", ori * 100 / (sampling.sample_points.size() - 1));
			for (int j = 1; j <= ori * 20 / sampling.sample_points.size(); j++)
				std::cout << "▉";

			Eigen::Vector3d vectorAfter;
			if (!PrepareOrientationSliceData(sampling, ori, now_last_node, save_ori, rotMatrix, vectorAfter)) {
				continue;
			}

			double temp_time_1 = 0;
			double temp_time_2 = 0;
			Layer_Graph layer_graph = BuildAdditiveLayerGraphWithSurfaceMesh(
				vectorAfter,
				height_of_beam_search,
				tree_entries[now_last_node].continue_id,
				the_nozzle,
				temp_time_1,
				temp_time_2);

			//Layer_Graph layer_graph = BuildAdditiveLayerGraph(
			//	vectorAfter,
			//	height_of_beam_search,
			//	tree_entries[now_last_node].continue_id,
			//	the_nozzle,
			//	sum_time_4,
			//	sum_time_2);

			//////////////////////////////////////
			vector<Eigen::MatrixXd> fragile_V;
			for (int i = 0; i < is_fragile_V_2[ori].size(); i++)
				if (is_fragile_V_2[ori][i] == true)
					fragile_V.push_back(V_2[i]);

			all_solutions_of_selected_layers.clear();
			all_solutions_of_selected_layers_contain.clear();
			all_cut_layers.clear(); all_cut_layers_dependency_layer.clear();
			pathes_include_S.clear(); pathes_include_sample_points.clear(); paths_include_covering_points.clear();
			//std::cout << "cc" << endl;


			bool flag_continue = false;	//是否继续从当前节点向下搜索
			bool jud_admit = true;
			if (!ProcessOrientationSearch(
				layer_graph,
				judge_S_be_searched,
				judge_covering_points_be_searched,
				W1,
				sum_time,
				flag_continue,
				jud_admit)) {
				continue;
			}

			final_pathes_include_S.insert(final_pathes_include_S.end(), pathes_include_S.begin(), pathes_include_S.end());
			final_pathes_include_sample_points.insert(final_pathes_include_sample_points.end(), pathes_include_sample_points.begin(), pathes_include_sample_points.end());
			final_pathes_include_covering_points.insert(final_pathes_include_covering_points.end(), paths_include_covering_points.begin(), paths_include_covering_points.end());

			for (int i = 0; i < all_solutions_of_selected_layers.size(); i++) {	//将当前方向生成的所有解加入树结构中，作为下一层搜索的节点
				index_node++;
				BeamTreeEntry entry;
				entry.layers = all_solutions_of_selected_layers[i];
				entry.contain = all_solutions_of_selected_layers_contain[i];
				entry.cut_layers = all_cut_layers[i];
				entry.cut_layers_dependency_layer = all_cut_layers_dependency_layer[i];
				entry.fragile_v = fragile_V;
				entry.parent_id = now_last_node;
				entry.pre_queue_index = cont_number_of_queue;
				candidate_nodes.push_back(index_node);
				vector<Eigen::Vector3d> temp_vec;
				temp_vec.push_back(vectorAfter);
				save_ori.push_back(temp_vec);
				entry.judge_continue = flag_continue;
				if (tree_entries[now_last_node].judge_continue == true)
					entry.continue_id = tree_entries[now_last_node].continue_id + 1;
				else
					entry.continue_id = 0;
				entry.error = false;
				tree_entries.push_back(entry);
			}
		}



		if (tree_entries[now_last_node].judge_continue == false)
			last_step_nodes.pop();
		cont_number_of_queue++;
		if (last_step_nodes.size() == 0 || tree_entries[now_last_node].judge_continue == true) {
			cont_number_of_queue = 0;
			height_of_beam_search++;
			int cont_w = 0;  //应为0开始
			cout << endl << "Decomposing and sorting......" << endl;
			/////////sort_candidate_nodes//////////
			start_time = clock();
			vector<all_value> all_calculated_value;
			vector<all_value> pure_value;
			all_calculated_value.resize(candidate_nodes.size());
			pure_value.resize(candidate_nodes.size());
			for (int i = 0; i < static_cast<int>(candidate_nodes.size()); i++) {
				EvaluateCandidateNode(
					i,
					candidate_nodes,
					all_calculated_value,
					tree_entries,
					save_ori,
					height_of_beam_search,
					cont_number_of_queue,
					file_name,
					now_last_node,
					the_nozzle);
			}

			sort_candidate_nodes(candidate_nodes, tree_entries, final_pathes_include_S, all_calculated_value, final_pathes_include_covering_points, height_of_beam_search, save_ori, pure_value, tree_entries[candidate_nodes[0]].continue_id);
			end_time = clock();

			ofile_cont << final_pathes_include_S[candidate_nodes[0]].size() << endl;
			double temp_record_size = 0;
			for (int j = 0; j < tree_entries[candidate_nodes[0]].layers.size(); j++) {
				for (int k = 0; k < tree_entries[candidate_nodes[0]].layers[j].size() - 1; k++) {
					temp_record_size += distance3d(tree_entries[candidate_nodes[0]].layers[j][k], tree_entries[candidate_nodes[0]].layers[j][k + 1]);
				}
			}
			ofile_cont_size << temp_record_size << endl;
			//////////////////////////////////////
			if (tree_entries[now_last_node].judge_continue == true) {
				W1 = 1;
			}

			//cout << "bbb";

			bool jud_continue_last_node = false;
			//	for (int i = 0; i < candidate_nodes.size(); i++)
			//		cout << save_ori[candidate_nodes[i]][0].x() << " " << save_ori[candidate_nodes[i]][0].y() << " " << save_ori[candidate_nodes[i]][0].z() << " " << endl;
			start_time_7 = clock();
			while (candidate_nodes.size() != 0 && cont_w < W1 && cont_w < candidate_nodes.size()) {
				LogSelectedCandidateMetrics(
					cont_w,
					candidate_nodes,
					now_last_node,
					save_ori,
					tree_entries,
					last_step_nodes,
					all_calculated_value,
					pure_value,
					Sum_candidate_blocks,
					Sum_connected_components,
					cont_extra_additive_orientation,
					evaluation_value);
				//decompose the model for every node//
				int index_of_pre_node = tree_entries[candidate_nodes[cont_w]].pre_queue_index;
				vector<vector<cv::Point3d>> all_cut_layers;
				vector<int> flag_cut_layers_is_hole;
				vector<int> all_cut_layers_dependency_layer;
				BuildCutLayersForCandidate(
					cont_w,
					candidate_nodes,
					tree_entries,
					all_cut_layers,
					all_cut_layers_dependency_layer,
					flag_cut_layers_is_hole);

				cout << "***** number_of_S: " << final_pathes_include_S[candidate_nodes[cont_w]].size() << endl;
				cout << "***** all_cut_layers.size(): " << all_cut_layers.size() << endl;

				std::vector<TRiangle> current_remove_triangles;
				Slicer_2 current_slicer;
				//cout << "aaa";
				/*for (int i = 0; i < Tree_nodes_continue_id.size(); i++)
					cout << "*"<<Tree_nodes_continue_id[i] << endl;*/

				CutMesh(tree_entries[candidate_nodes[cont_w]].layers, tree_entries[candidate_nodes[cont_w]].contain,
					all_cut_layers, save_ori[candidate_nodes[cont_w]][0],
					height_of_beam_search,
					cont_number_of_queue,
					index_of_pre_node,
					all_cut_layers_dependency_layer,
					jud_outer_beam_search_terminate,
					current_remove_triangles,
					current_slicer,
					tree_entries[candidate_nodes[cont_w]].judge_continue, tree_entries[now_last_node].judge_continue,
					tree_entries[now_last_node].pre_queue_index,
					tree_entries[candidate_nodes[cont_w]].error,
					candidate_nodes[cont_w],
					tree_entries[candidate_nodes[cont_w]].continue_id,
					flag_cut_layers_is_hole);


				//if (Tree_nodes_judge_continue[now_last_node] == true) {
				//	cout << "aaa";
				//	Tree_nodes.pop_back();
				//	Tree_nodes_cut_layers.pop_back();
				//	Tree_nodes_num_of_cut_layers_dependency_layer.pop_back();
				//	Tree_nodes_fragile_V.pop_back();
				//	index_node--;
				//	candidate_nodes.pop_back();
				//	save_ori.pop_back();
				//	Tree_nodes_judge_continue.pop_back();
				//	Tree_nodes_error.pop_back();
				//	height_of_beam_search--;
				//	Tree_nodes_judge_continue[now_last_node] = false;
				//	//Tree_nodes_continue_id[now_last_node]++;
				//}
				cout << "&&&&" << height_of_beam_search << endl;
				if (tree_entries[candidate_nodes[0]].judge_continue == true) {
					cout << "AAAA!" << endl;
					jud_continue_last_node = true;

				}
				if (tree_entries[now_last_node].judge_continue == true) { //Tree_nodes_judge_continue[now_last_node] == true  //应该改为Tree_nodes_judge_continue[candidate_nodes[0]] == true?
					cout << "BBBB!" << candidate_nodes[0] << endl;
					if (tree_entries[candidate_nodes[0]].judge_continue == false)
						tree_entries[now_last_node].judge_continue = false;
					height_of_beam_search--;
					tree_entries.pop_back();
					index_node--;
					candidate_nodes.pop_back();
					//Tree_nodes_judge_continue[now_last_node] = false;
					tree_entries[now_last_node].continue_id = tree_entries[candidate_nodes[0]].continue_id;
					save_ori[now_last_node].push_back(save_ori[candidate_nodes[0]][0]);
					save_ori.pop_back();

				}


				if (height_of_beam_search <= 150) {
					//cout << save_ori[candidate_nodes[cont_w]].x() <<" "<< save_ori[candidate_nodes[cont_w]].y() << " " << save_ori[candidate_nodes[cont_w]].z() << endl;
					subtractive_accessibility_decomposition(current_remove_triangles, height_of_beam_search, cont_number_of_queue, cutting_tool, current_slicer);
					//subtractive_accessibility_decomposition(current_remove_triangles, height_of_beam_search, index_of_pre_node, cutting_tool, current_slicer);
				}

				subtractive_remove_output(current_remove_triangles, current_slicer, height_of_beam_search, cont_number_of_queue);

				cont_number_of_queue++;
				//////////////////////////////////////

				cont_w++;
			}
			end_time_7 = clock();
			std::cout << "&&&time&&& Evaluation " << double(end_time - start_time + end_time_7 - start_time_7) / CLOCKS_PER_SEC << std::endl;
			std::cout << "&&&time&&& Update subtractive dependency graph: " << sum_time_3 << std::endl;
			std::cout << "&&&time&&& Slicing: " << sum_time_4 << std::endl;
			std::cout << "&&&time&&& Build addictive dependency graph: " << sum_time_2 << std::endl;
			std::cout << "&&&time&&& Co graph merging: " << sum_time << std::endl;
			total_time_1 += sum_time_3;
			total_time_2 += sum_time_2 + sum_time_4 + sum_time_5;
			total_time_3 += sum_time;
			total_time_4 += double(end_time - start_time + end_time_7 - start_time_7) / CLOCKS_PER_SEC;
			sum_time = 0;
			sum_time_2 = 0;
			sum_time_3 = 0;
			sum_time_4 = 0;
			sum_time_5 = 0;
			//std::cout << "&&&time&&& a step of beam search: " << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << std::endl;
			/*if (height_of_beam_search == 5)
				W1 = 20;*/
				/*if (height_of_beam_search >= 2 && Tree_nodes_continue_id[candidate_nodes[0]] >= 0) {
					end_time_total = clock();
					std::cout << "&&&time&&& total of outer-beam search: " << double(end_time_total - start_time_total) / CLOCKS_PER_SEC << "s &&&&&&&" << std::endl << endl;
					system("pause");
				} */

				//system("pause");
				//save_ori.clear();
			cont_number_of_queue = 0;
			if (jud_outer_beam_search_terminate == true) {
				/*cout << "$$$$$$$$$$$$$$$$$$$$$$$$$ Data $$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
				cout << "#S: " << all_voronoi_cells.size() << endl;
				cout << "#P_i: " << num_inaccessible_points << endl;
				cout << "#R_p: " << float(num_inaccessible_points) / float(all_voronoi_cells.size()) << endl;
				cout << "#B_C !!: " << Sum_candidate_blocks << endl;
				cout << "#I: " << height_of_beam_search << endl;
				cout << "#B: " << height_of_beam_search << endl;
				cout << "#O: " << cont_extra_additive_orientation + height_of_beam_search << endl;
				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$ Total time $$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
				cout << "Build subtractive graph: " << time_build_subtractive_graph << endl;
				cout << "Update subtractive dependency graph: " << total_time_1 << endl;
				cout << "Build addictive dependency graph: " << total_time_2 << endl;
				cout << "Co graph merging: " << total_time_3 << endl;
				cout << "Evaluation: " << total_time_4 << endl;
				cout << "Total: " << time_build_subtractive_graph + total_time_1 + total_time_2 + total_time_3 + total_time_4 << endl;
				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;*/

				double global_score = 1 / (height_of_beam_search);
				fstream f;
				//追加写入,在原来基础上加了ios::app  faith
				std::string txt_name = file_name + ".txt";
				f.open(txt_name.c_str(), ios::out | ios::app);
				//输入你想写入的内容 
				f << global_score << endl;
				f.close();

				break;
			}
			if (height_of_beam_search >= 10) {
				fstream f;
				//追加写入,在原来基础上加了ios::app  faith
				std::string txt_name = file_name + ".txt";
				f.open(txt_name.c_str(), ios::out | ios::app);
				//输入你想写入的内容 
				f << 0 << endl;
				f.close();
				//break;
			}
			if (jud_continue_last_node == true)
				cout << endl << "--------------------- height of outer-beam search: " << height_of_beam_search << " ---------------------";
			else
				cout << endl << "--------------------- height of outer-beam search: " << height_of_beam_search + 1 << " ---------------------";
			candidate_nodes.clear();
		}
	}

	vector<int> last_available_block;
	vector<int> exist_points;
	vector<vector<int>> is_point_exist_in_block;
	vector<vector<int>> is_point_available_in_block;
	vector<int> is_available(height_of_beam_search + 1, 1);


	last_available_block.resize(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		last_available_block[i] = height_of_beam_search;
	}

	for (int i = 0; i < V.rows(); i++) {
		is_point_exist_in_block.push_back(is_available);
	}

	for (int i = 0; i < V.rows(); i++)
		is_point_available_in_block.push_back(is_available);

	cout << "TTYTYTYTYTY";

	//subtractive_accessibility_decomposition_global(height_of_beam_search, cutting_tool);
}

double HybridManufacturing::UpdateSubtractiveDependencyGraph(
	const vector<vector<area_S>>& ori_all_the_area_S,
	vector<vector<bool>>& is_fragile_V_2,
	int now_last_node,
	int height_of_beam_search,
	int cont_number_of_queue,
	const vector<BeamTreeEntry>& tree_entries,
	vector<bool>& judge_S_be_searched,
	vector<bool>& judge_covering_points_be_searched)
{
	clock_t start_time_5 = clock();
	vector<int> S_in_block;
	vector<int> sample_points_in_block;	//记录不存在于当前块中的顶点索引
	vector<int> covering_points_in_block;	//记录被碰撞单元索引
	judge_S_be_searched.clear();
	judge_covering_points_be_searched.clear();
	judge_S_be_searched.resize(all_the_area_S.size());
	judge_covering_points_be_searched.resize(all_the_covering_points.size());
	for (int i = 0; i < all_the_area_S.size(); i++)
		judge_S_be_searched[i] = false;
	for (int i = 0; i < all_the_covering_points.size(); i++)
		judge_covering_points_be_searched[i] = false;
	S_in_block.clear();
	sample_points_in_block.clear();
	covering_points_in_block.clear();
	for (int i = 0; i < V.rows(); i++) {	//枚举原始输入模型网格V的所有顶点，判断该顶点是否仍存在于当前块中
		bool jud_still_exist = false;
		/*if (flag_sample_point_used[i])
			continue;*/
		for (int j = 0; j < Katana::Instance().vertices.size(); j++) {
			if (abs(V(i, 0) - Katana::Instance().vertices[j].x) <= 0.001 && abs(V(i, 1) - Katana::Instance().vertices[j].y) <= 0.001 && abs(V(i, 2) - Katana::Instance().vertices[j].z) <= 0.001) {
				jud_still_exist = true;	//标记该顶点仍存在
				//flag_sample_point_used[i] = true;
				break;
			}
		}

		if (jud_still_exist == false) {	//如果该顶点不存在于当前块中
			for (int k = 0; k < is_fragile_V_2.size(); k++)
				is_fragile_V_2[k][i] = false;	//is_fragile_V_2中对应该顶点的信息全部置为false
			sample_points_in_block.push_back(i);
			for (int k = 0; k < map_S_and_vertex.size(); k++) {
				if (map_S_and_vertex[k] == i)
					S_in_block.push_back(k);	//将该顶点对应的不可达索引加入S_in_block
			}
			for (int k = 0; k < map_covering_points_and_vertex.size(); k++) {
				if (map_covering_points_and_vertex[k] == i)
					covering_points_in_block.push_back(k);	//将该顶点对应的被碰撞单元索引加入covering_points_in_block
			}
		}
	}

	/*while (pre_tree_nodes[now_last_node] != -1) {
		for (int i = 0; i < final_pathes_include_S[now_last_node].size(); i++) {
			S_in_block.push_back(final_pathes_include_S[now_last_node][i]);
		}
		for (int i = 0; i < final_pathes_include_sample_points[now_last_node].size(); i++) {
			sample_points_in_block.push_back(final_pathes_include_sample_points[now_last_node][i]);
		}
		now_last_node = pre_tree_nodes[now_last_node];
	}
	now_last_node = last_step_nodes.front();*/

	for (int i = 0; i < S_in_block.size(); i++)
		judge_S_be_searched[S_in_block[i]] = true;

	ori_num_points_of_ori_in_all_the_area_S.assign(
		ori_all_the_area_S.size(),
		vector<int>(sampling_subtractive.sample_points.size(), 0)); //为每个不可达点在ori_num_points_of_ori_in_all_the_area_S分配一个向量
	//每个不可达点记录在各个采样方向上对应的点数

	for (size_t point_id = 0; point_id < ori_all_the_area_S.size(); ++point_id) {
		const auto& orientations_for_point = ori_all_the_area_S[point_id];
		auto& orientation_point_count = ori_num_points_of_ori_in_all_the_area_S[point_id];
		for (const auto& area : orientations_for_point) {
			++orientation_point_count[area.oriId];
		}
	}

	for (int missing_vertex_id : sample_points_in_block) {
		const int covering_point_vertex_id = map_covering_points_and_vertex_inv[missing_vertex_id];
		const auto& covering_points_for_vertex = all_the_covering_points[covering_point_vertex_id];
		for (const auto& covering_point : covering_points_for_vertex) {
			ori_num_points_of_ori_in_all_the_area_S[covering_point.pointId][covering_point.oriId]--;
		}
	}

	for (int covering_point_id : covering_points_in_block) {
		judge_covering_points_be_searched[covering_point_id] = true;
	}

	if (tree_entries[now_last_node].judge_continue == true) {
		Katana::Instance().stl.loadStl((file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(tree_entries[now_last_node].continue_id) + "_subblock.stl").c_str());
		Katana::Instance().temp_vertices.resize(Katana::Instance().vertices.size());
		Katana::Instance().temp_vertices = Katana::Instance().vertices;
	}
	clock_t end_time_5 = clock();
	return double(end_time_5 - start_time_5) / CLOCKS_PER_SEC;
}

bool HybridManufacturing::PrepareOrientationSliceData(
	const SAMPLE_ON_BALL& sampling,
	int ori,
	int now_last_node,
	const vector<vector<Eigen::Vector3d>>& save_ori,
	Eigen::Matrix3d& rot_matrix,
	Eigen::Vector3d& vector_after)
{
	Eigen::Vector3d vector_continue(sampling.sample_points[ori]);  //与祖先Node的ori相同时容易出现问题，暂时先避免用相同ori
	bool jud_continue = false;

	for (int j = 0; j < save_ori[now_last_node].size(); j++)
		if (vector_continue == save_ori[now_last_node][j]) {
			jud_continue = true;
			break;
		}

	if (jud_continue == true && ori != 0) {
		return false;
	}

	std::vector<Eigen::MatrixXd> temp_V;	//记录当前块的顶点信息，用于旋转
	temp_V.resize(Katana::Instance().temp_vertices.size());

	V_2.resize(V.rows());	//V_2记录输入模型网格V的顶点信息，用于旋转
	for (int i = 0; i < Katana::Instance().temp_vertices.size(); i++) {
		temp_V[i].resize(3, 1);
		temp_V[i](0, 0) = Katana::Instance().temp_vertices[i].x;
		temp_V[i](1, 0) = Katana::Instance().temp_vertices[i].y;
		temp_V[i](2, 0) = Katana::Instance().temp_vertices[i].z;
	}
	for (int i = 0; i < V.rows(); i++) {
		V_2[i].resize(3, 1);
		V_2[i](0, 0) = V.row(i).x();
		V_2[i](1, 0) = V.row(i).y();
		V_2[i](2, 0) = V.row(i).z();
	}
	Eigen::Vector3d vector_before(0, 0, 1);
	vector_after = Eigen::Vector3d(sampling.sample_points[ori]);
	rot_matrix = Eigen::Quaterniond::FromTwoVectors(vector_before, vector_after).toRotationMatrix();
	for (int i = 0; i < Katana::Instance().vertices.size(); i++)
		temp_V[i] = rot_matrix.inverse() * temp_V[i];	//将temp_V里的顶点坐标旋转，使得增材打印方向与z轴平行
	for (int i = 0; i < V.rows(); i++)
		V_2[i] = rot_matrix.inverse() * V_2[i];	//将V_2里的顶点坐标旋转，使得增材打印方向与z轴平行
	for (int i = 0; i < Katana::Instance().vertices.size(); i++) {
		Katana::Instance().vertices[i].x = temp_V[i](0, 0);
		Katana::Instance().vertices[i].y = temp_V[i](1, 0);
		Katana::Instance().vertices[i].z = temp_V[i](2, 0);	//将temp_V里的顶点坐标更新回Katana的vertices
	}

	for (int i = 0; i < Katana::Instance().triangles.size(); i++) {	//将每个三角形的顶点按z值从小到大排序，便于后续分层处理
		if (Katana::Instance().triangles[i].vertices[0]->z > Katana::Instance().triangles[i].vertices[1]->z) std::swap(Katana::Instance().triangles[i].vertices[0], Katana::Instance().triangles[i].vertices[1]);
		if (Katana::Instance().triangles[i].vertices[0]->z > Katana::Instance().triangles[i].vertices[2]->z) std::swap(Katana::Instance().triangles[i].vertices[0], Katana::Instance().triangles[i].vertices[2]);
		if (Katana::Instance().triangles[i].vertices[1]->z > Katana::Instance().triangles[i].vertices[2]->z) std::swap(Katana::Instance().triangles[i].vertices[1], Katana::Instance().triangles[i].vertices[2]);
	}

	current_node_mesh_rotated.clear();
	current_node_mesh_rotated = current_node_mesh;

	for (auto it : current_node_mesh_rotated.vertices())
	{
		auto& vertex = current_node_mesh_rotated.point(it);
		auto temp_vec = Eigen::Vector3d(vertex.x(), vertex.y(), vertex.z());
		temp_vec = rot_matrix.inverse() * temp_vec;
		vertex = Point_3(temp_vec.x(), temp_vec.y(), temp_vec.z());
	}

	CGAL::IO::write_polygon_mesh("model\\temp_current_node_mesh.obj", current_node_mesh);
	CGAL::IO::write_polygon_mesh("model\\temp_rotated.obj", current_node_mesh_rotated);

	return true;
}

// 保持函数签名与声明一致，避免不必要的声明差异。
Layer_Graph HybridManufacturing::BuildAdditiveLayerGraph(
	const Eigen::Vector3d& vector_after,
	int height_of_beam_search,
	int continue_node_id,
	const nozzle& the_nozzle,
	double& slicing_time,
	double& graph_time)
{
	clock_t start_time_6 = clock();
	vector<vector<vector<Vertex>>> all_slice_points; //(i,j,k) 第i层第j个loop的第k个点的坐标
	vector<vector<vector<Vertex>>> all_slice_points_contain;
	Katana::Instance().slicer.buildLayers();


	bool segments_manifold = Katana::Instance().slicer.buildSegments();
	if (!segments_manifold) {
		std::cout << "[HybridManufacturing::BuildAdditiveLayerGraph] !segments_manifold" << std::endl;
		//Visualize_layer_segments(Katana::Instance().layers);

		CGAL::Polygon_mesh_slicer<SurfaceMesh, Kernel> slicer_cgal(current_node_mesh_rotated);
		std::vector<Polylines> layer_polylines;
		layer_polylines.resize(Katana::Instance().layers.size());
		for (auto& lay : Katana::Instance().layers) {
			slicer_cgal(Kernel::Plane_3(0, 0, -1, lay.z), std::back_inserter(layer_polylines[&lay - &Katana::Instance().layers[0]]));
		}

		//Visualize_layer_polylines(layer_polylines);
	}

	Katana::Instance().gcode.write(all_slice_points, all_slice_points_contain);	//katana分层并得到边界轮廓

	vector<vector<vector<Eigen::Vector3d>>> vector_after_all_slice_points;	//记录每个切片轮廓上每个点的坐标
	vector<vector<vector<Eigen::Vector3d>>> vector_after_all_slice_points_contain;
	vector_after_all_slice_points.resize(all_slice_points.size());
	vector_after_all_slice_points_contain.resize(all_slice_points_contain.size());
	for (int i = 0; i < all_slice_points.size(); i++) {
		vector_after_all_slice_points[i].resize(all_slice_points[i].size());
		for (int j = 0; j < all_slice_points[i].size(); j++) {
			vector_after_all_slice_points[i][j].resize(all_slice_points[i][j].size());
			for (int k = 0; k < all_slice_points[i][j].size(); k++) {
				vector_after_all_slice_points[i][j][k] = Eigen::Vector3d(all_slice_points[i][j][k].x, all_slice_points[i][j][k].y, all_slice_points[i][j][k].z);
			}
		}
	}
	for (int i = 0; i < all_slice_points_contain.size(); i++) {
		vector_after_all_slice_points_contain[i].resize(all_slice_points_contain[i].size());
		for (int j = 0; j < all_slice_points_contain[i].size(); j++) {
			vector_after_all_slice_points_contain[i][j].resize(all_slice_points_contain[i][j].size());
			for (int k = 0; k < all_slice_points_contain[i][j].size(); k++) {
				vector_after_all_slice_points_contain[i][j][k] = Eigen::Vector3d(all_slice_points_contain[i][j][k].x, all_slice_points_contain[i][j][k].y, all_slice_points_contain[i][j][k].z);
			}
		}
	}

	clock_t end_time_6 = clock();
	slicing_time += double(end_time_6 - start_time_6) / CLOCKS_PER_SEC;
	/////////////////////////////////////
	//std::cout << "aa" << Katana::Instance().vertices.size()<<" " << all_slice_points.size() << endl;
	//generate additive dependency graph//
	std::vector<Data> data;
	data.resize(1);
	data[0].ReadData(vector_after_all_slice_points, vector_after_all_slice_points_contain);
	Layer_Graph layer_graph(data[0]);
	clock_t start_time_4 = clock();
	layer_graph.GetTrianglesForLayers(all_slice_points, Katana::Instance().map_segment_triangles, Katana::Instance().vertices, vector_after, height_of_beam_search, continue_node_id);	//将切片轮廓映射到三角形集，建立每层三角形集合等中间信息？
	layer_graph.GenerateDependencyEdges();	//生成增材分层依赖图的边
	layer_graph.CollisionDetectionForAdditiveManufacturing(the_nozzle);	//增材的碰撞检测，标记增材的不可达点等信息
	clock_t end_time_4 = clock();
	graph_time += double(end_time_4 - start_time_4) / CLOCKS_PER_SEC;

	return layer_graph;
}

bool HybridManufacturing::ProcessOrientationSearch(
	Layer_Graph& layer_graph,
	const vector<bool>& judge_S_be_searched,
	const vector<bool>& judge_covering_points_be_searched,
	int& W1,
	double& sum_time,
	bool& flag_continue,
	bool& jud_admit)
{
	bool previous_is_continue = false;	//记录上一个节点是否为continue节点
	clock_t start_time_3 = clock();

	if (open_change_orientation == true)
		W1 = 1;
	DFS_search(layer_graph, flag_continue, previous_is_continue, judge_S_be_searched, judge_covering_points_be_searched, jud_admit);	//对当前方向进行增材分层依赖图的深度优先搜索，生成这个分层方向的解

	/*if (height_of_beam_search==1&& Tree_nodes_continue_id[now_last_node] >= 2)
		flag_continue = false;*/

	clock_t end_time_3 = clock();
	if (jud_admit == false)
		return false;
	sum_time += double(end_time_3 - start_time_3) / CLOCKS_PER_SEC;

	//cout << "()()()(" << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << endl;
	return true;
}

void HybridManufacturing::DFS_search(Layer_Graph layer_graph, bool& flag_continue, bool previous_is_continue, vector<bool> judge_S_be_searched, vector<bool> judge_covering_points_be_searched, bool& jud_admit)
{
	int W2 = 1;  //2
	struct DfsTreeEntry {
		int layer = -1;
		vector<int> area_s;
		vector<int> sample_points;
		vector<int> covering_points;
		int parent = -1;
	};
	vector<DfsTreeEntry> tree_entries;
	vector<int> candidate_nodes;
	queue<int> last_step_nodes;
	vector<int> S_in_block;
	vector<int> terminate_nodes;
	DfsTreeEntry root_entry;
	tree_entries.push_back(root_entry);
	int index_node = 0;
	last_step_nodes.push(index_node);
	//vector<vector<area_S>> ori_all_the_area_S = all_the_area_S;

	//std::cout << "AA" << endl;
	//find boundary point of every layer//
	vector<double> left_point_of_layer, right_point_of_layer, top_point_of_layer, bottom_point_of_layer;
	left_point_of_layer.resize(layer_graph.total_node_num); right_point_of_layer.resize(layer_graph.total_node_num);
	top_point_of_layer.resize(layer_graph.total_node_num); bottom_point_of_layer.resize(layer_graph.total_node_num);

	vector<vector<Point_2>> vec_points_2d(layer_graph.total_node_num);
	vector<Polygon_2> vec_polygon(layer_graph.total_node_num);
	/*Point_2 pp(V(i, 0), V(i, 1));
	if (polygon.bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE)
		jud_still_exist = true;*/

	for (int i = 0; i < layer_graph.total_node_num; i++) {
		left_point_of_layer[i] = MAX_I;
		right_point_of_layer[i] = -MAX_I;
		top_point_of_layer[i] = -MAX_I;
		bottom_point_of_layer[i] = MAX_I;
	}
	int temp_num = 0;
	for (int i = 0; i < layer_graph.data.slice_points.size(); i++) {
		for (int j = 0; j < layer_graph.data.slice_points[i].size(); j++) {
			for (int k = 0; k < layer_graph.data.slice_points[i][j].size(); k++) {
				/*if (layer_graph.data.slice_points[i][j][k].x < left_point_of_layer[temp_num])
					left_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].x;
				if (layer_graph.data.slice_points[i][j][k].x > right_point_of_layer[temp_num])
					right_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].x;*/
					/*if (layer_graph.data.slice_points[i][j][k].y < bottom_point_of_layer[temp_num])
						bottom_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].y;
					if (layer_graph.data.slice_points[i][j][k].y > top_point_of_layer[temp_num])
						top_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].y;*/
				vec_points_2d[temp_num].push_back(Point_2(layer_graph.data.slice_points[i][j][k].x(), layer_graph.data.slice_points[i][j][k].y()));
			}
			vec_polygon[temp_num] = Polygon_2(vec_points_2d[temp_num].begin(), vec_points_2d[temp_num].end());
			//std::cout << "!!@@##" << std::endl;
			if (!vec_polygon[temp_num].is_simple()) {

				vec_polygon[temp_num] = Polygon_2(vec_points_2d[temp_num].begin(), vec_points_2d[temp_num].end() - 1);

				if (!vec_polygon[temp_num].is_simple()) {
					std::cout << "polygon is not simple HybridManufacturing::DFS_search cutcut!" << std::endl;
					std::string poly_file_name = file_name + "not_simple" + std::to_string(temp_num) + ".txt";
					std::ofstream poly_file(poly_file_name);
					poly_file << "Polygon points count: " << vec_points_2d[temp_num].size() << std::endl;
					for (const auto& point : vec_points_2d[temp_num]) {
						poly_file << point.x() << " " << point.y() << std::endl;
					}
					poly_file.close();
				}
			}
			temp_num++;
		}
	}


	//////////////////////////////////////
	//std::cout << "BB" << endl;
	//find all the sample points in every layer//
	vector<vector<int>> sample_point_in_layer;	//存储每个layer中包含的sample point的索引
	sample_point_in_layer.resize(layer_graph.data.total_node_num);
	vector<bool> judge_sample_point_be_searched;
	judge_sample_point_be_searched.resize(V_2.size());
	for (int i = 0; i < V_2.size(); i++)
		judge_sample_point_be_searched[i] = false;
	for (int i = 0; i < layer_graph.data.total_node_num; i++) {
		for (int j = 0; j < V_2.size(); j++) {
			pair<int, int> index_slice_layer = layer_graph.data.index[i];
			Point_2 pp(V_2[j](0, 0), V_2[j](1, 0));
			if (judge_sample_point_be_searched[j] == false
				&& vec_polygon[i].bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE // bounded_sideh很耗时?
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[j](2, 0) > 0
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[j](2, 0) <= dh) {
				sample_point_in_layer[i].push_back(j);	//V_2的第j个点在第i个layer中
				judge_sample_point_be_searched[j] = true;	//标记该点需要被搜索
			}
		}
	}
	////////////////////////////////////////

	//std::cout << "CC" << endl;
	//find area S and covering points in the node//
	vector<vector<int>> temp_all_S_in_the_block, temp_all_covering_points_in_the_block;
	temp_all_S_in_the_block.clear(); temp_all_covering_points_in_the_block.clear();
	temp_all_S_in_the_block.resize(layer_graph.data.total_node_num);
	temp_all_covering_points_in_the_block.resize(layer_graph.data.total_node_num);
	for (int i = 0; i < layer_graph.data.total_node_num; i++) {
		pair<int, int> index_slice_layer = layer_graph.data.index[i];
		for (int j = 0; j < ori_num_points_of_ori_in_all_the_area_S.size(); j++) {
			Point_2 pp(V_2[map_S_and_vertex[j]](0, 0), V_2[map_S_and_vertex[j]](1, 0));
			if (judge_S_be_searched[j] == false
				&& vec_polygon[i].bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[map_S_and_vertex[j]](2, 0) > 0
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[map_S_and_vertex[j]](2, 0) <= dh) {
				temp_all_S_in_the_block[i].push_back(j);
				judge_S_be_searched[j] = true;	//标记该area S需要被搜索
			}
		}
		for (int j = 0; j < ori_all_the_covering_points.size(); j++) {
			Point_2 pp(V_2[map_covering_points_and_vertex[j]](0, 0), V_2[map_covering_points_and_vertex[j]](1, 0));
			if (judge_covering_points_be_searched[j] == false
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[map_covering_points_and_vertex[j]](2, 0) > 0
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[map_covering_points_and_vertex[j]](2, 0) <= dh
				&& vec_polygon[i].bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE) {
				temp_all_covering_points_in_the_block[i].push_back(j);
				judge_covering_points_be_searched[j] = true;	//标记该covering point需要被搜索
			}
		}
	}
	////////////////////////////////////////
	//std::cout << "DD" << endl;
	while (last_step_nodes.size() != 0 || candidate_nodes.size() != 0) {
		flag_continue = false;
		int now_last_node = last_step_nodes.front();
		bool jud_terminate = true;

		//update//
		while (tree_entries[now_last_node].parent != -1) {
			layer_graph.UpdateDegree_2(tree_entries[now_last_node].layer, -1);
			layer_graph.node_visited[tree_entries[now_last_node].layer] = true;
			now_last_node = tree_entries[now_last_node].parent;
		}
		//////////
		now_last_node = last_step_nodes.front();

		for (int i = 0; i < layer_graph.total_node_num; i++) {
			int v = i;
			bool flag_self_support = true;
			if (tree_entries.size() != 0 && tree_entries[now_last_node].layer == v) continue;
			bool jud_continue_2 = false;
			/*for (int j = 0; j < candidate_nodes.size(); j++)
				if (Tree_nodes[candidate_nodes[j]] == v) {
					jud_continue_2 = true;
					break;
				}*/
			if (jud_continue_2 == true) continue;
			if (layer_graph.node_visited[v]) continue;
			if (layer_graph.out_degree[v] != 0) continue;   //out-degree == 0
			//////////collision dependency edges///////////////
			if (tree_entries.size() != 0 && layer_graph.IsDepend_collision(tree_entries[now_last_node].layer, v))
				continue;
			///////////////////////////////////////////////////
			if (layer_graph.is_the_layer_self_suppot[v] == false && now_last_node != 0) {  //self-support constraint
				flag_self_support = false;
			}

			vector<int> all_S_in_the_block, all_covering_points_in_the_block;
			all_S_in_the_block = temp_all_S_in_the_block[i];
			all_covering_points_in_the_block = temp_all_covering_points_in_the_block[i];

			//try merge//
			bool jud_merge = true;
			for (int j = 0; j < all_S_in_the_block.size(); j++) {
				//check whether all S in the block is accessible//
				bool jud_merge_2 = false;
				for (int k = 0; k < ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]].size(); k++) {
					if (ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]][k] < 0)
						cout << "error!!";
					if (ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]][k] == 0) {
						jud_merge_2 = true;
						break;
					}
				}
				if (jud_merge_2 == false)
					jud_merge = false;
				//////////////////////////////////////////////////
			}
			if (jud_merge == true && flag_self_support == false && previous_is_continue == false) {
				if (open_change_orientation == true)
					flag_continue = true;
			}
			if (jud_merge == true && flag_self_support == true) {
				jud_terminate = false;
				index_node++;
				DfsTreeEntry entry;
				entry.layer = v;
				entry.area_s = all_S_in_the_block;
				entry.sample_points = sample_point_in_layer[v];
				entry.covering_points = all_covering_points_in_the_block;
				entry.parent = now_last_node;
				tree_entries.push_back(entry);
				candidate_nodes.push_back(index_node);

				break; //DFS
			}
			////////////
		}
		if (jud_terminate == true) {
			terminate_nodes.push_back(now_last_node);
			break;
		}

		//update//
		while (tree_entries[now_last_node].parent != -1) {
			layer_graph.UpdateDegree_2(tree_entries[now_last_node].layer, 1);
			layer_graph.node_visited[tree_entries[now_last_node].layer] = false;
			now_last_node = tree_entries[now_last_node].parent;
		}
		//////////

		//sort_candidate_nodes(candidate_nodes, Tree_nodes_for_S);
		last_step_nodes.pop();
		//cout << candidate_nodes.size() << endl;
		if (last_step_nodes.size() == 0) {
			int cont_w = 0;
			while (candidate_nodes.size() != 0 && cont_w < W2 && cont_w < candidate_nodes.size()) {
				last_step_nodes.push(candidate_nodes[cont_w]);
				cont_w++;
			}
			candidate_nodes.clear();
		}
	}

	//every path (blcok) as a candidate node of outer-beam search (save slice layers)//
	//final_pathes save the index of slice layer (layer_graph.total_node_num)
	pathes_include_S.clear();
	pathes_include_sample_points.clear();
	paths_include_covering_points.clear();
	vector<vector<int>> final_pathes;
	vector<int> final_pathes_height;
	for (int i = 0; i < terminate_nodes.size(); i++) {
		int height = 0;
		vector<int> current_path, current_path_include_S, current_path_include_sample_points, current_path_include_covering_points;
		current_path.clear();
		current_path_include_S.clear();
		current_path_include_sample_points.clear();
		current_path_include_covering_points.clear();
		int current_node = terminate_nodes[i];
		current_path.push_back(tree_entries[terminate_nodes[i]].layer);
		for (int j = 0; j < tree_entries[terminate_nodes[i]].area_s.size(); j++)
			current_path_include_S.push_back(tree_entries[terminate_nodes[i]].area_s[j]);
		for (int j = 0; j < tree_entries[terminate_nodes[i]].sample_points.size(); j++)
			current_path_include_sample_points.push_back(tree_entries[terminate_nodes[i]].sample_points[j]);
		for (int j = 0; j < tree_entries[terminate_nodes[i]].covering_points.size(); j++)
			current_path_include_covering_points.push_back(tree_entries[terminate_nodes[i]].covering_points[j]);
		while (tree_entries[current_node].parent > 0) {
			current_node = tree_entries[current_node].parent;
			current_path.push_back(tree_entries[current_node].layer);
			height++;
			for (int j = 0; j < tree_entries[current_node].area_s.size(); j++) {
				current_path_include_S.push_back(tree_entries[current_node].area_s[j]);
			}

			for (int j = 0; j < tree_entries[current_node].sample_points.size(); j++)
				current_path_include_sample_points.push_back(tree_entries[current_node].sample_points[j]);
			for (int j = 0; j < tree_entries[current_node].covering_points.size(); j++)
				current_path_include_covering_points.push_back(tree_entries[current_node].covering_points[j]);
		}
		std::reverse(current_path.begin(), current_path.end());
		std::reverse(current_path_include_S.begin(), current_path_include_S.end());
		std::reverse(current_path_include_sample_points.begin(), current_path_include_sample_points.end());
		std::reverse(current_path_include_covering_points.begin(), current_path_include_covering_points.end());
		final_pathes.push_back(current_path);
		final_pathes_height.push_back(height);
		pathes_include_S.push_back(current_path_include_S);
		pathes_include_sample_points.push_back(current_path_include_sample_points);
		paths_include_covering_points.push_back(current_path_include_covering_points);
	}

	//sort final_pathes and select W2 pathes//
	for (int i = 0; i < final_pathes.size(); i++)
		for (int j = i + 1; j < final_pathes.size(); j++) {
			if (final_pathes_height[i] < final_pathes_height[j]) {
				swap(final_pathes[i], final_pathes[j]);
				swap(final_pathes_height[i], final_pathes_height[j]);
			}
		}
	while (final_pathes.size() > W2) {
		auto itr = final_pathes.begin() + W2;
		final_pathes.erase(itr);
		auto itr2 = pathes_include_S.begin() + W2;
		pathes_include_S.erase(itr2);
		auto itr3 = pathes_include_sample_points.begin() + W2;
		pathes_include_sample_points.erase(itr3);
		auto itr4 = paths_include_covering_points.begin() + W2;
		paths_include_covering_points.erase(itr4);
	}
	////////////////////

	all_cut_layers = FindAllCutLayers(layer_graph, final_pathes, all_cut_layers_dependency_layer, jud_admit);

	vector<vector<vector<cv::Point3d>>> all_selected_points;	//存储所有被选择的路径上的切片轮廓的各个点的坐标
	vector<vector<vector<cv::Point3d>>> all_selected_points_contain;
	for (int i = 0; i < final_pathes.size(); i++) {
		vector<vector<cv::Point3d>> temp_vec_1;
		all_selected_points.push_back(temp_vec_1);
		all_selected_points_contain.push_back(temp_vec_1);
		for (int j = 0; j < final_pathes[i].size(); j++) {
			vector<cv::Point3d> temp_vec_2;
			all_selected_points[i].push_back(temp_vec_2);
			all_selected_points_contain[i].push_back(temp_vec_2);
			pair<int, int> index_slice_point = layer_graph.data.index[final_pathes[i][j]];
			for (int k = 0; k < layer_graph.data.slice_points[index_slice_point.first][index_slice_point.second].size(); k++) {
				cv::Point3d current_point(layer_graph.data.slice_points[index_slice_point.first][index_slice_point.second][k].x(),
					layer_graph.data.slice_points[index_slice_point.first][index_slice_point.second][k].y(),
					layer_graph.data.z_value[index_slice_point.first][index_slice_point.second][k]);
				all_selected_points[i][j].push_back(current_point);
			}
			for (int k = 0; k < layer_graph.data.slice_points_contain[index_slice_point.first][index_slice_point.second].size(); k++) {
				cv::Point3d current_point(layer_graph.data.slice_points_contain[index_slice_point.first][index_slice_point.second][k].x(),
					layer_graph.data.slice_points_contain[index_slice_point.first][index_slice_point.second][k].y(),
					layer_graph.data.z_value[index_slice_point.first][index_slice_point.second][0]);
				all_selected_points_contain[i][j].push_back(current_point);
			}
		}
	}


	//////////visualization all triangles cross by layers//////////
	//vector<vector<vector<Vertex*>>> show_triangles(final_pathes[0].size());
	//for (int i = 0; i < 1; i++) {
	//	for (int j = 0; j < final_pathes[0].size(); j++) {
	//		show_triangles[j].resize(layer_graph.all_triangles_of_layers[final_pathes[i][j]].size());
	//		for (int k = 0; k < layer_graph.all_triangles_of_layers[final_pathes[i][j]].size(); k++) {
	//			Vertex* v1 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[0];
	//			Vertex* v2 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[1];
	//			Vertex* v3 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[2];
	//			//double ans = (v2->x - v1->x) * (v2->y - v3->y) - (v2->y - v1->y) * (v2->x - v3->x);
	//			//if (ans > 0)	//is clockwise
	//			//	swap(v2, v3);
	//			show_triangles[j][k].push_back(v1);
	//			show_triangles[j][k].push_back(v2);
	//			show_triangles[j][k].push_back(v3);
	//		}
	//	}
	//}
	//std::ofstream dstream(".\\vis\\show_triangles.stl");
	//if (!dstream.is_open()) {
	//	std::cout << "can not open " << std::endl;
	//}
	//dstream << "solid STL generated by MeshLab" << std::endl;
	//for (int j = 0; j < final_pathes[0].size(); j++) {
	//	for (int k = 0; k < show_triangles[j].size(); k++) {
	//		dstream << "  facet normal " << "0 0 0" << std::endl;
	//		dstream << "    outer loop" << std::endl;
	//		for (int l = 0; l < 3; l++) {
	//			dstream << "      vertex  " << show_triangles[j][k][l]->x << " " << show_triangles[j][k][l]->y << " " << show_triangles[j][k][l]->z << std::endl;
	//		}
	//		dstream << "    endloop" << std::endl;
	//		dstream << "  endfacet" << std::endl;
	//	}
	//}
	//dstream << "endsolid vcg" << std::endl;
	//dstream.close();

	all_solutions_of_selected_layers = all_selected_points;
	all_solutions_of_selected_layers_contain = all_selected_points_contain;
	return;
	//////////////////////////////////////////////////////////////////////////////////
}

void HybridManufacturing::sort_candidate_nodes(
	vector<int>& candidate_nodes,
	const vector<BeamTreeEntry>& tree_entries,
	vector<vector<int>> final_pathes_include_S,
	OrientationScores& all_calculated_value,
	vector<vector<int>> final_pathes_include_covering_points,
	int height_of_beam_search,
	vector<vector<Eigen::Vector3d>> save_ori,
	OrientationScores& pure_value,
	int id_continue)
{
	double W_less_area_S = 0, W_more_slice_layers = 0.6, W_covering_points = 0.4, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;

	//double W_less_area_S = 0, W_more_slice_layers = 0.1, W_covering_points = 0.1, W_less_clipping_plane = 0.4, W_fragile = 0, W_larger_base = 0, W_orientation = 0.4, W_projected = 0;
	//double W_less_area_S = 0, W_more_slice_layers = 0, W_covering_points = 0, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 1;

	//if(height_of_beam_search >= 12)  //test for figure of paper
	//	W_less_area_S = 0.1, W_more_slice_layers = 1, W_covering_points = 0.1, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;
	/*if(height_of_beam_search == 2)
		W_less_area_S = 0.9, W_more_slice_layers = 0.1, W_covering_points = 0, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;*/

		//random test
		//int rand_id = rand() % all_calculated_value.size();
		/*int rand_id = 5;
		swap(candidate_nodes[0], candidate_nodes[rand_id]);
		swap(all_calculated_value[0], all_calculated_value[rand_id]);
		return;*/

	for (int i = 0; i < all_calculated_value.size(); i++) {
		if (all_calculated_value[i].number_of_remaining_face < terminate_threshold_of_number_of_faces && save_ori[candidate_nodes[i]][0].x() == 0 && save_ori[candidate_nodes[i]][0].y() == 0) {
			swap(candidate_nodes[0], candidate_nodes[i]);
			swap(all_calculated_value[0], all_calculated_value[i]);
			cout << "Remaining face: " << all_calculated_value[i].number_of_remaining_face << endl;
			return;
		}
	}
	for (int i = 0; i < all_calculated_value.size(); i++) {
		if (height_of_beam_search == 4 && save_ori[candidate_nodes[i]][0].x() == 0 && save_ori[candidate_nodes[i]][0].y() == 0) {
			swap(candidate_nodes[0], candidate_nodes[i]);
			swap(all_calculated_value[0], all_calculated_value[i]);
			cout << "Remaining face: " << all_calculated_value[i].number_of_remaining_face << endl;
			return;
		}
	}

	/*if(height_of_beam_search == 2)
		W_less_area_S = 0, W_more_slice_layers = 0.1, W_covering_points = 0.1, W_less_clipping_plane = 0.8, W_fragile = 0, W_larger_base = 0;*/

		//normalization - local//
		/*int max_area_S = -MAX_I, min_area_S = MAX_I;
		int max_slice_layers = -MAX_I, min_slice_layers = MAX_I;
		int max_clipping_plane = -MAX_I, min_clipping_plane = MAX_I;
		double max_large_base = -MAX_I, min_large_base = MAX_I;
		for (int i = 0; i < all_calculated_value.size(); i++) {
			max_large_base = max(max_large_base, all_calculated_value[i].large_base);
			min_large_base = min(min_large_base, all_calculated_value[i].large_base);

			max_area_S = max(max_area_S, int(final_pathes_include_S[candidate_nodes[i]].size()));
			min_area_S = min(min_area_S, int(final_pathes_include_S[candidate_nodes[i]].size()));

			max_slice_layers = max(max_slice_layers, int(tree_entries[candidate_nodes[i]].layers.size()));
			min_slice_layers = min(min_slice_layers, int(tree_entries[candidate_nodes[i]].layers.size()));

			max_clipping_plane = max(max_clipping_plane, int(tree_entries[candidate_nodes[i]].cut_layers.size()));
		}


		for (int i = 0; i < all_calculated_value.size(); i++) {
			if (max_large_base - min_large_base != 0)
				all_calculated_value[i].large_base = (all_calculated_value[i].large_base - min_large_base) / (max_large_base - min_large_base);
			else
				all_calculated_value[i].large_base = 0;

			if (max_area_S - min_area_S != 0)
				all_calculated_value[i].value_of_area_S = double(final_pathes_include_S[candidate_nodes[i]].size() - min_area_S) / double(max_area_S - min_area_S);
			else
				all_calculated_value[i].value_of_area_S = 0;

			if (max_slice_layers - min_slice_layers != 0)
			all_calculated_value[i].value_of_more_slice_layers = double(tree_entries[candidate_nodes[i]].layers.size() - min_slice_layers) / double(max_slice_layers - min_slice_layers);
			else
				all_calculated_value[i].value_of_more_slice_layers = 0;

			if (max_clipping_plane - min_clipping_plane != 0)
			all_calculated_value[i].value_of_less_clipping_plane = 1 - double(tree_entries[candidate_nodes[i]].cut_layers.size() - min_clipping_plane) / double(max_clipping_plane - min_clipping_plane);
			else
				all_calculated_value[i].value_of_less_clipping_plane = 0;
		}*/
		/////////////////


		//Calculate the accumulated value//
	int* sum_area_S = new int[all_calculated_value.size()];
	double* sum_slice_layers = new double[all_calculated_value.size()];
	int* sum_clipping_plane = new int[all_calculated_value.size()];
	double* sum_larger_base = new double[all_calculated_value.size()];
	double* sum_covering_points = new double[all_calculated_value.size()];
	double* sum_fragile = new double[all_calculated_value.size()];
	double* sum_orientation = new double[all_calculated_value.size()];
	double* sum_projected = new double[all_calculated_value.size()];
	for (int i = 0; i < all_calculated_value.size(); i++) {
		sum_area_S[i] = 0;
		sum_slice_layers[i] = 0;
		sum_clipping_plane[i] = 0;
		sum_larger_base[i] = 0;
		sum_covering_points[i] = 0;
		sum_orientation[i] = 0;
		sum_fragile[i] = all_calculated_value[i].value_of_fragile;
		//cout << "((((((" << sum_fragile[i] << endl;
		sum_projected[i] = all_calculated_value[i].value_of_projected;
		int index_node = candidate_nodes[i];
		const int parent_id = tree_entries[index_node].parent_id;
		if (parent_id >= 0) {
			pure_value[i].value_of_orientation = save_ori[index_node][0].dot(save_ori[parent_id][0]);
		}
		else {
			pure_value[i].value_of_orientation = 0;
		}
		pure_value[i].value_of_fragile = all_calculated_value[i].value_of_fragile;
		pure_value[i].value_of_projected = all_calculated_value[i].value_of_projected;
		while (tree_entries[index_node].parent_id != -1) {
			sum_area_S[i] += int(final_pathes_include_S[index_node].size());
			for (int j = 0; j < tree_entries[index_node].layers.size(); j++) {
				for (int k = 0; k < tree_entries[index_node].layers[j].size() - 1; k++) {
					sum_slice_layers[i] += distance3d(tree_entries[index_node].layers[j][k], tree_entries[index_node].layers[j][k + 1]);
				}
			}
			sum_clipping_plane[i] += int(tree_entries[index_node].cut_layers.size());
			sum_larger_base[i] += tree_entries[index_node].larger_base;

			//sum_covering_points[i] += int(final_pathes_include_covering_points[index_node].size());
			for (int j = 0; j < final_pathes_include_covering_points[index_node].size(); j++)
				sum_covering_points[i] += int(all_the_covering_points[final_pathes_include_covering_points[index_node][j]].size());
			index_node = tree_entries[index_node].parent_id;
		}
		index_node = candidate_nodes[i];
		while (tree_entries[index_node].parent_id > 0) {
			sum_orientation[i] += save_ori[index_node][0].dot(save_ori[tree_entries[index_node].parent_id][0]);
			index_node = tree_entries[index_node].parent_id;
		}
	}
	////////////////////////////////////


	//normalization - accumulated//
	int max_area_S = -MAX_I, min_area_S = MAX_I;
	double max_slice_layers = -MAX_I, min_slice_layers = MAX_I;
	int max_clipping_plane = -MAX_I, min_clipping_plane = MAX_I;
	double max_large_base = -MAX_I, min_large_base = MAX_I;
	double max_covering_points = -MAX_I, min_covering_points = MAX_I;
	double max_fragile = -MAX_I, min_fragile = MAX_I;
	double max_orientation = -MAX_I, min_orientation = MAX_I;
	double max_projected = -MAX_I, min_projected = MAX_I;
	for (int i = 0; i < all_calculated_value.size(); i++) {
		max_large_base = max(max_large_base, sum_larger_base[i]);
		min_large_base = min(min_large_base, sum_larger_base[i]);

		max_area_S = max(max_area_S, sum_area_S[i]);
		min_area_S = min(min_area_S, sum_area_S[i]);

		max_slice_layers = max(max_slice_layers, sum_slice_layers[i]);
		min_slice_layers = min(min_slice_layers, sum_slice_layers[i]);

		max_clipping_plane = max(max_clipping_plane, sum_clipping_plane[i]);
		min_clipping_plane = min(min_clipping_plane, sum_clipping_plane[i]);

		max_covering_points = max(max_covering_points, sum_covering_points[i]);
		min_covering_points = min(min_covering_points, sum_covering_points[i]);

		max_fragile = max(max_fragile, sum_fragile[i]);
		min_fragile = min(min_fragile, sum_fragile[i]);

		max_orientation = max(max_orientation, sum_orientation[i]);
		min_orientation = min(min_orientation, sum_orientation[i]);

		max_projected = max(max_projected, sum_projected[i]);
		min_projected = min(min_projected, sum_projected[i]);
	}

	for (int i = 0; i < all_calculated_value.size(); i++) {
		if (max_large_base - min_large_base != 0)
			all_calculated_value[i].large_base = (sum_larger_base[i] - min_large_base) / (max_large_base - min_large_base);
		else
			all_calculated_value[i].large_base = 1;  //全部暂时改为1，之后还是需要-1用于标识

		if (max_area_S - min_area_S != 0)
			all_calculated_value[i].value_of_area_S = 1 - double(sum_area_S[i] - min_area_S) / double(max_area_S - min_area_S);
		else
			all_calculated_value[i].value_of_area_S = 1;

		if (max_slice_layers - min_slice_layers != 0)
			all_calculated_value[i].value_of_more_slice_layers = double(sum_slice_layers[i] - min_slice_layers) / double(max_slice_layers - min_slice_layers);
		else
			all_calculated_value[i].value_of_more_slice_layers = 1;

		if (max_clipping_plane - min_clipping_plane != 0)
			all_calculated_value[i].value_of_less_clipping_plane = 1 - double(sum_clipping_plane[i] - min_clipping_plane) / double(max_clipping_plane - min_clipping_plane);
		else
			all_calculated_value[i].value_of_less_clipping_plane = 1;

		if (max_covering_points - min_covering_points != 0)
			all_calculated_value[i].value_of_covering_points = double(sum_covering_points[i] - min_covering_points) / double(max_covering_points - min_covering_points);
		else
			all_calculated_value[i].value_of_covering_points = 1;

		if (max_fragile - min_fragile != 0)
			all_calculated_value[i].value_of_fragile = 1 - double(sum_fragile[i] - min_fragile) / double(max_fragile - min_fragile);
		else
			all_calculated_value[i].value_of_fragile = 1;

		if (max_orientation - min_orientation != 0)
			all_calculated_value[i].value_of_orientation = double(sum_orientation[i] - min_orientation) / double(max_orientation - min_orientation);
		else
			all_calculated_value[i].value_of_orientation = 1;

		/*if (height_of_beam_search == 3) {
			if (max_projected - min_projected != 0)
				all_calculated_value[i].value_of_projected = double(sum_projected[i] - min_projected) / double(max_projected - min_projected);
			else
				all_calculated_value[i].value_of_projected = -1;
		}*/
		//else {
		if (max_projected - min_projected != 0)
			all_calculated_value[i].value_of_projected = 1 - double(sum_projected[i] - min_projected) / double(max_projected - min_projected);
		else
			all_calculated_value[i].value_of_projected = 1;
		//}
	}
	/////////////////

	//sort
	for (int i = 0; i < candidate_nodes.size(); i++) {
		for (int j = i + 1; j < candidate_nodes.size(); j++) {
			auto calc_score = [&](int index) {
				const auto& value = all_calculated_value[index];
				return value.value_of_area_S * W_less_area_S
					+ value.value_of_more_slice_layers * W_more_slice_layers
					+ value.value_of_less_clipping_plane * W_less_clipping_plane
					+ value.large_base * W_larger_base
					+ value.value_of_covering_points * W_covering_points
					+ value.value_of_fragile * W_fragile
					+ value.value_of_orientation * W_orientation
					+ value.value_of_projected * W_projected
					+ value.value_of_sub_patches * 0.5;
				};

			if (calc_score(i) < calc_score(j)) {
				swap(candidate_nodes[i], candidate_nodes[j]);
				swap(all_calculated_value[i], all_calculated_value[j]);
				swap(pure_value[i], pure_value[j]);
			}
		}
	}


	//if (height_of_beam_search == 6) {
		//cout << endl << "&&" << all_calculated_value[0].number_of_remaining_face << endl;
		//cout << endl << "&&" << all_calculated_value[1].number_of_remaining_face << endl;
		//}

	//先删除
	cout << "#####################" << height_of_beam_search << " " << id_continue << "#####################" << endl;
	//if (height_of_beam_search == 2) {  //julia_vase
	//	int cont = 0;
	//	while (cont < 2) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	// 
	//if (height_of_beam_search == 3) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//else if (height_of_beam_search == 5) {
	//	int cont = 0;
	//	while (cont <1) {  //6
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 2 && id_continue == 0) {
	//	int cont = 0;
	//	while (cont < 1) {  //6
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}

	//if (height_of_beam_search == 4) {
	//	int cont = 0;
	//	while (cont < 5) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 5) {
	//	int cont = 0;
	//	while (cont < 8) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 6) {
	//	int cont = 0;
	//	while (cont < 13) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 7) {
	//	int cont = 0;
	//	while (cont < 0) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 13) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 14) {
	//	int cont = 0;
	//	while (cont < 2) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 18) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
}

void HybridManufacturing::sort_candidate_nodes(vector<int>& candidate_nodes, vector<vector<int>> Tree_nodes_for_S)
{
	/*for (int i = 0; i < candidate_nodes.size(); i++) {
		for (int j = i + 1; j < candidate_nodes.size(); j++) {
			if (Tree_nodes_for_S[candidate_nodes[i]].size() < Tree_nodes_for_S[candidate_nodes[j]].size()) {
				swap(candidate_nodes[i], candidate_nodes[j]);
			}
		}
	}*/
}

void HybridManufacturing::subtractive_remove_output(const vector<TRiangle>& need_detect_triangle, const Slicer_2& current_slicer, int height_of_beam_search, int cont_number_of_queue)
{
	std::string filename = ".\\vis\\block_patch-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".obj";
	std::ofstream file(filename);
	if (!file)
	{
		std::cout << "subtractive_remove_output !file" << std::endl;
	}

	for (auto& p : current_slicer.positions) {
		file << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
	}
	for (auto t : need_detect_triangle) {
		for (int i = 0; i < 3; ++i) t[i]++;
		file << "f " << t[0] << " " << t[1] << " " << t[2] << "\n";
	}
}

bool HybridManufacturing::CheckToolCollisionWithCell(
	const Eigen::Vector3d& center_point,
	const std::vector<Eigen::Vector3d>& target_cell_vertices,
	double max_z_target,
	const cutter& nozzle,
	double z_threshold_divisor,
	double xy_tolerance) const
{
	(void)xy_tolerance;

	const double z_threshold = nozzle.cylinder_r / z_threshold_divisor;
	const double height_diff = max_z_target - center_point.z();

	// 快速排除：目标在刀尖球下方
	if (height_diff <= z_threshold) {
		return false;
	}

	// 快速确认：目标超出刀具最大高度
	if (height_diff > nozzle.total_height) {
		return true;
	}

	// 使用第一个顶点进行粗筛
	if (!target_cell_vertices.empty()) {
		const double dx = target_cell_vertices[0](0, 0) - center_point.x();
		const double dy = target_cell_vertices[0](1, 0) - center_point.y();
		const double dist_xy_sq = dx * dx + dy * dy;

		if (height_diff > nozzle.cylinder_height_threshold) {
			if (dist_xy_sq > nozzle.carriage_check_radius_sq) {
				return false;
			}
		}
		else {
			if (dist_xy_sq > nozzle.cylinder_check_radius_sq) {
				return false;
			}
		}
	}

	// 精确检测每个边界顶点
	for (const auto& vertex : target_cell_vertices) {
		const double vx = vertex(0, 0);
		const double vy = vertex(1, 0);
		const double vz = vertex(2, 0);
		const double diff_z = vz - center_point.z();

		if (diff_z <= z_threshold) {
			continue;
		}

		const double dx = vx - center_point.x();
		const double dy = vy - center_point.y();
		const double dist_xy_sq = dx * dx + dy * dy;

		if (diff_z <= nozzle.cylinder_height_threshold) { //nozzle.cylinder_r + nozzle.cylinder_height
			if (dist_xy_sq < nozzle.cylinder_r_sq) {
				return true;
			}
		}
		else if (diff_z <= nozzle.total_height) { //nozzle.cylinder_r + nozzle.cylinder_height + nozzle.carriage_height
			if (dist_xy_sq < nozzle.carriage_r_sq) {
				return true;
			}
		}
	}

	return false;
}
