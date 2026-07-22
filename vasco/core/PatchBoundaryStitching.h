#pragma once

#include <cstddef>
#include <vector>

#include "vasco/Slicer.h"

namespace vasco::patch_boundary_stitching {

struct StitchOptions {
	double model_relative_tolerance = 1e-6;
	double edge_relative_tolerance = 1e-3;
	double absolute_min_tolerance = 1e-9;
	double absolute_max_tolerance = 1e-3;
	double normal_cosine_threshold = 0.95;
	double boundary_direction_cosine_threshold = 0.98;
	int max_iterations = 8;
};

struct StitchStats {
	double model_diagonal = 0.0;
	double model_tolerance = 0.0;
	std::size_t initial_boundary_edge_count = 0;
	std::size_t initial_t_junction_count = 0;
	std::size_t projected_vertex_count = 0;
	std::size_t split_edge_count = 0;
	std::size_t split_face_count = 0;
	std::size_t added_face_count = 0;
	std::size_t rejected_candidate_count = 0;
	std::size_t remaining_t_junction_count = 0;
	std::size_t final_boundary_edge_count = 0;
	std::size_t iteration_count = 0;
};

// Stitches cross-patch boundary T-junctions by projecting newer boundary
// vertices onto older boundary edges and splitting the older incident faces.
// Split triangles inherit the original face source id and preserve its winding.
bool StitchPatchBoundariesWithTolerance(
	vasco::Slicer& mesh,
	std::vector<int>& face_source_patch_ids,
	const StitchOptions& options,
	StitchStats& stats);

} // namespace vasco::patch_boundary_stitching
