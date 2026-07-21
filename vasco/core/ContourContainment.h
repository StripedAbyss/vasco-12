#pragma once

#include <vector>

#include "surface_mesh_slice_data.h"

namespace vasco::contour_containment
{
	// One printable material region represented by an even-depth outer contour
	// and its direct odd-depth hole contours.
	struct ContourComponent {
		int outer_contour_id = -1;
		std::vector<int> hole_contour_ids;
	};

	// Builds a parent-only containment tree for contours in one slice.
	// Each result entry is the smallest enclosing contour, or -1 for a root.
	std::vector<int> BuildContourParentIds(
		const Polylines& contours,
		double coordinate_tolerance = 1e-8);

	// Derives the direct-child list for every node from the stored parent IDs.
	std::vector<std::vector<int>> BuildDirectChildren(
		const std::vector<int>& parent_ids);

	// Counts parent links from a contour to its root; returns -1 for invalid data.
	int ComputeContourDepth(
		const std::vector<int>& parent_ids,
		int contour_id);

	// Groups every even-depth contour with all of its direct odd-depth holes.
	std::vector<ContourComponent> BuildMaterialComponents(
		const std::vector<int>& parent_ids);
}
