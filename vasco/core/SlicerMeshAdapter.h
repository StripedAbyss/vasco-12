#pragma once

#include <vector>

#include "surface_mesh_slice_data.h"
#include "vasco/Slicer.h"

namespace vasco::slicer_mesh_adapter {

void RemeshCutPatchBeforeSaving(
	vasco::Slicer& all_slicer,
	std::vector<int>& id_contact_faces,
	const std::vector<double>& cut_plane_z_values,
	double target_edge_length);

} // namespace vasco::slicer_mesh_adapter
