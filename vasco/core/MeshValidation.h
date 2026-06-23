#pragma once

#include <string>

#include "surface_mesh_slice_data.h"

namespace vasco::mesh_validation {

bool LoadAndCheckCurrentNodeMesh(const std::string& current_file_name, SurfaceMesh& current_node_mesh);

} // namespace vasco::mesh_validation
