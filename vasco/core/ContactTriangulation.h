#pragma once

#include <map>
#include <vector>

#include "surface_mesh_slice_data.h"
#include "vasco/Slicer.h"
#include "vasco/core/Types.h"

namespace vasco::contact_triangulation {

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

QuantizedPointKey MakeQuantizedKey(const Point_3& p, double scale = 1e9);

std::vector<vasco::core::Tri3> TriangulateContactFacesWithCDT(
	const std::vector<int>& outer_ring,
	const std::vector<std::vector<int>>& hole_rings,
	vasco::Slicer& slicer,
	std::vector<int>& id_contact_faces,
	std::map<QuantizedPointKey, int>& point_index_map);

} // namespace vasco::contact_triangulation
