#include "vasco/core/ContactTriangulation.h"

#include <algorithm>
#include <cmath>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

namespace vasco::contact_triangulation {
namespace {

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
	ContactFaceBase(
		typename Fb::Vertex_handle v0,
		typename Fb::Vertex_handle v1,
		typename Fb::Vertex_handle v2,
		typename Fb::Face_handle n0,
		typename Fb::Face_handle n1,
		typename Fb::Face_handle n2)
		: Fb(v0, v1, v2, n0, n1, n2) {}

private:
	bool in_domain_ = false;
};

using ContactFb = ContactFaceBase<ContactTriangulationTraits>;
using ContactTds = CGAL::Triangulation_data_structure_2<ContactVb, ContactFb>;
using ContactCdt = CGAL::Constrained_Delaunay_triangulation_2<ContactTriangulationTraits, ContactTds>;

double SignedArea2D(const Point_3& a, const Point_3& b, const Point_3& c)
{
	return CGAL::to_double((b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()));
}

int GetOrCreateSlicerPositionIndex(
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

void AddConstraintLoop(ContactCdt& cdt, const std::vector<int>& ring, const vasco::Slicer& slicer)
{
	if (ring.size() < 3) return;

	std::vector<ContactCdt::Vertex_handle> handles;
	handles.reserve(ring.size());
	for (int idx : ring) {
		const auto& pos = slicer.positions[idx];
		handles.push_back(cdt.insert(ContactPoint2(pos[0], pos[1])));
	}
	for (std::size_t i = 0; i < handles.size(); ++i) {
		cdt.insert_constraint(handles[i], handles[(i + 1) % handles.size()]);
	}
}

double ComputeAverageBoundaryEdgeLength(const std::vector<int>& ring, const vasco::Slicer& slicer)
{
	if (ring.size() < 2) {
		return 1.0;
	}

	double sum = 0.0;
	for (std::size_t i = 0; i < ring.size(); ++i) {
		const auto& a = slicer.positions[ring[i]];
		const auto& b = slicer.positions[ring[(i + 1) % ring.size()]];
		double dx = a[0] - b[0];
		double dy = a[1] - b[1];
		sum += std::sqrt(dx * dx + dy * dy);
	}
	return std::max(sum / static_cast<double>(ring.size()), 1e-6);
}

} // namespace

QuantizedPointKey MakeQuantizedKey(const Point_3& p, double scale)
{
	return {
		static_cast<long long>(std::llround(CGAL::to_double(p.x()) * scale)),
		static_cast<long long>(std::llround(CGAL::to_double(p.y()) * scale)),
		static_cast<long long>(std::llround(CGAL::to_double(p.z()) * scale))
	};
}

std::vector<vasco::core::Tri3> TriangulateContactFacesWithCDT(
	const std::vector<int>& outer_ring,
	const std::vector<std::vector<int>>& hole_rings,
	vasco::Slicer& slicer,
	std::vector<int>& id_contact_faces,
	std::map<QuantizedPointKey, int>& point_index_map)
{
	std::vector<vasco::core::Tri3> generated;
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
	using Criteria = CGAL::Delaunay_mesh_size_criteria_2<ContactCdt>;
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

		vasco::core::Tri3 tri{};
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

} // namespace vasco::contact_triangulation
