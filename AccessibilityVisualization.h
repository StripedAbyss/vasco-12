#pragma once

#include <Eigen/Dense>

#include <string>
#include <unordered_map>
#include <vector>

#include "helpers.h"
#include "katana/datastructures.h"
#include "surface_mesh_slice_data.h"
#include "vasco/Slicer.h"
#include "vasco/core/Voronoi.h"

class Layer_Graph;

namespace accessibility_visualization
{
	// Writes the initial set of points that are inaccessible in every subtractive direction.
	void WriteInitialSubtractiveAllDirectionDebugVisualization(
		const std::string& vis_dir,
		const std::string& file_token,
		double marker_radius,
		const Eigen::MatrixXd& original_vertices,
		const std::unordered_map<int, int>& map_s_to_vertex);

	// Writes one tool/blocker OBJ per orientation for the highest active inaccessible point.
	void WriteHighestZInaccessiblePointAllOrientationToolCollisionObj(
		const std::string& output_file,
		const Eigen::MatrixXd& original_vertices,
		const std::vector<vasco::VoronoiCell>& voronoi_cells,
		const std::vector<Eigen::Vector3d>& orientation_samples,
		const std::unordered_map<int, int>& map_s_to_vertex,
		cutter tool,
		int& selected_s_id,
		int& selected_cell_id,
		const std::vector<bool>* active_cell_mask = nullptr,
		const std::vector<bool>* searched_s_flags = nullptr);

	// Writes the per-point collision-reason table and its tool-collision visualizations.
	void WriteSubtractiveAllDirectionReasonDiagnostics(
		const std::string& vis_dir,
		const std::string& file_token,
		const Eigen::MatrixXd& original_vertices,
		const std::vector<vasco::VoronoiCell>& voronoi_cells,
		const std::vector<Eigen::Vector3d>& orientation_samples,
		const std::unordered_map<int, int>& map_s_to_vertex,
		cutter tool);

	// Writes the unresolved subtractive-S and covering-point marker meshes for one beam node.
	void WriteSubtractiveAccessibilityDebugVisualizations(
		const std::string& vis_dir,
		const std::string& node_tag,
		double marker_radius,
		const Eigen::MatrixXd& original_vertices,
		const std::vector<bool>& searched_s_flags,
		const std::vector<bool>& searched_covering_point_flags,
		const std::unordered_map<int, int>& map_s_to_vertex,
		const std::unordered_map<int, int>& map_covering_point_to_vertex);

	// Marks original vertices that are still present in the current cut-node mesh.
	std::vector<bool> BuildActiveOriginalVertexMask(
		const Eigen::MatrixXd& original_vertices,
		const std::vector<Vertex>& current_vertices,
		double eps = 0.001);

	// Writes marker meshes for additive layer nodes targeted by collision-dependency edges.
	void WriteAdditiveAccessibilityDebugVisualization(
		const std::string& vis_dir,
		const Layer_Graph& layer_graph,
		const std::string& node_tag,
		int orientation_id,
		double marker_radius);

	// Writes the out-degree-zero DFS root candidates with the minimum/maximum-z layers as references.
	void WriteAdditiveRootLayerSelfSupportDebugVisualization(
		const std::string& vis_dir,
		const Layer_Graph& layer_graph,
		const std::string& node_tag,
		int orientation_id,
		const Eigen::Vector3d& orientation,
		const SurfaceMesh& mesh);

	// Shows the final graph-cut labels as one Polyscope mesh and tool-direction arrow per label.
	void ShowMergedPatchGraphCutLabels(
		const std::string& vis_dir,
		const std::string& output_prefix,
		int height_of_beam_search,
		int orientation_count,
		const vasco::Slicer& merged_patch,
		const std::vector<int>& result_labels,
		const std::vector<Eigen::Vector3d>& orientation_samples);
}
