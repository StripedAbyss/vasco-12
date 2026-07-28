#include "AccessibilityVisualization.h"

#include "layer_graph.h"
#include "visual.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <direct.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>

namespace accessibility_visualization
{
	namespace
	{
		// Precomputes the cutter dimensions used by the collision diagnostics.
		void PrepareToolForCollision(cutter& tool)
		{
			tool.cylinder_height_threshold = tool.cylinder_height + tool.ball_r;
			tool.carriage_check_radius_sq = (tool.carriage_r + 5.0) * (tool.carriage_r + 5.0);
			tool.cylinder_check_radius_sq = (tool.cylinder_r + 5.0) * (tool.cylinder_r + 5.0);
			tool.cylinder_r_sq = tool.cylinder_r * tool.cylinder_r;
			tool.carriage_r_sq = tool.carriage_r * tool.carriage_r;
			tool.total_height = tool.cylinder_height + tool.ball_r + tool.carriage_height;
		}

		// Computes the tool-axis center used by the production collision test.
		Eigen::Vector3d ComputeToolCenter(
			const Eigen::Vector3d& point,
			const Eigen::Vector3d& normal,
			double radius)
		{
			return point + radius * normal;
		}

		// Creates the immediate parent directory required by a debug output file.
		void EnsureParentDirectory(const std::string& file_path)
		{
			const std::size_t slash_pos = file_path.find_last_of("/\\");
			if (slash_pos != std::string::npos) {
				_mkdir(file_path.substr(0, slash_pos).c_str());
			}
		}

		// Joins one directory and file name using the Windows path separator.
		std::string JoinPath(const std::string& directory, const std::string& file_name)
		{
			if (directory.empty()) {
				return file_name;
			}
			const char last = directory.back();
			if (last == '\\' || last == '/') {
				return directory + file_name;
			}
			return directory + "\\" + file_name;
		}
	}

	void WriteDebugMarkersObj(
		const std::string& output_file,
		const std::vector<Eigen::Vector3d>& points,
		const std::array<double, 3>& color,
		double marker_radius)
	{
		EnsureParentDirectory(output_file);
		std::ofstream ofs(output_file);
		if (!ofs.is_open()) {
			std::cout << "[AccessibilityDebug] cannot open file for writing: " << output_file << std::endl;
			return;
		}

		ofs << std::setprecision(17);
		ofs << "# marker_count " << points.size() << "\n";
		const std::array<Eigen::Vector3d, 6> offsets = {
			Eigen::Vector3d(marker_radius, 0.0, 0.0),
			Eigen::Vector3d(-marker_radius, 0.0, 0.0),
			Eigen::Vector3d(0.0, marker_radius, 0.0),
			Eigen::Vector3d(0.0, -marker_radius, 0.0),
			Eigen::Vector3d(0.0, 0.0, marker_radius),
			Eigen::Vector3d(0.0, 0.0, -marker_radius)
		};

		int vertex_base = 1;
		for (const auto& point : points) {
			for (const auto& offset : offsets) {
				const Eigen::Vector3d v = point + offset;
				ofs << "v " << v.x() << " " << v.y() << " " << v.z() << " "
					<< color[0] << " " << color[1] << " " << color[2] << "\n";
			}

			const int xp = vertex_base;
			const int xn = vertex_base + 1;
			const int yp = vertex_base + 2;
			const int yn = vertex_base + 3;
			const int zp = vertex_base + 4;
			const int zn = vertex_base + 5;
			ofs << "f " << zp << " " << xp << " " << yp << "\n";
			ofs << "f " << zp << " " << yp << " " << xn << "\n";
			ofs << "f " << zp << " " << xn << " " << yn << "\n";
			ofs << "f " << zp << " " << yn << " " << xp << "\n";
			ofs << "f " << zn << " " << yp << " " << xp << "\n";
			ofs << "f " << zn << " " << xn << " " << yp << "\n";
			ofs << "f " << zn << " " << yn << " " << xn << "\n";
			ofs << "f " << zn << " " << xp << " " << yn << "\n";
			vertex_base += static_cast<int>(offsets.size());
		}

		std::cout << "[AccessibilityDebug] wrote " << points.size()
			<< " markers to " << output_file << std::endl;
	}

	bool GetLayerGraphNodeCentroid(const Layer_Graph& layer_graph, int node_id, Eigen::Vector3d& centroid)
	{
		const auto index_it = layer_graph.data.index.find(node_id);
		if (index_it == layer_graph.data.index.end()) {
			return false;
		}

		const int layer_id = index_it->second.first;
		const int contour_id = index_it->second.second;
		if (layer_id < 0 || layer_id >= static_cast<int>(layer_graph.data.slice_points.size())
			|| contour_id < 0 || contour_id >= static_cast<int>(layer_graph.data.slice_points[layer_id].size())
			|| contour_id >= static_cast<int>(layer_graph.data.z_value[layer_id].size())) {
			return false;
		}

		const auto& contour = layer_graph.data.slice_points[layer_id][contour_id];
		const auto& z_values = layer_graph.data.z_value[layer_id][contour_id];
		if (contour.empty() || z_values.empty()) {
			return false;
		}

		Eigen::Vector2d xy_sum(0.0, 0.0);
		for (const auto& point : contour) {
			xy_sum += point;
		}
		xy_sum /= static_cast<double>(contour.size());
		centroid = Eigen::Vector3d(xy_sum.x(), xy_sum.y(), z_values.front());
		return true;
	}

	std::vector<Eigen::Vector3d> CollectMappedVertexPoints(
		const Eigen::MatrixXd& original_vertices,
		const std::unordered_map<int, int>& point_to_vertex)
	{
		std::vector<Eigen::Vector3d> points;
		points.reserve(point_to_vertex.size());
		for (const auto& entry : point_to_vertex) {
			const int vertex_id = entry.second;
			if (vertex_id < 0 || vertex_id >= original_vertices.rows()) {
				continue;
			}

			points.emplace_back(
				original_vertices(vertex_id, 0),
				original_vertices(vertex_id, 1),
				original_vertices(vertex_id, 2));
		}
		return points;
	}

	void WriteInitialSubtractiveAllDirectionDebugVisualization(
		const std::string& vis_dir,
		const std::string& file_token,
		double marker_radius,
		const Eigen::MatrixXd& original_vertices,
		const std::unordered_map<int, int>& map_S_and_vertex)
	{
		const auto inaccessible_points = CollectMappedVertexPoints(original_vertices, map_S_and_vertex);
		WriteDebugMarkersObj(
			JoinPath(vis_dir, "access_debug_subtractive_all_direction_inaccessible_" + file_token + ".obj"),
			inaccessible_points,
			{ 0.95, 0.05, 0.85 },
			marker_radius);
		WriteDebugMarkersObj(
			JoinPath(vis_dir, "access_debug_subtractive_all_direction_inaccessible_centers_" + file_token + ".obj"),
			inaccessible_points,
			{ 0.95, 0.05, 0.85 },
			0.03);
	}

	struct ToolCollisionDebugInfo {
		std::string reason = "none";
		int vertex_index = -1;
		double height_diff = 0.0;
		double dist_xy_sq = 0.0;
		double z_threshold = 0.0;
	};

	bool CheckToolCollisionWithCellForDebug(
		const Eigen::Vector3d& center_point,
		const std::vector<Eigen::Vector3d>& target_cell_vertices,
		double max_z_target,
		const cutter& tool,
		ToolCollisionDebugInfo& debug_info,
		double z_threshold_divisor = 30.0)
	{
		debug_info = ToolCollisionDebugInfo{};
		debug_info.z_threshold = tool.cylinder_r / z_threshold_divisor;
		debug_info.height_diff = max_z_target - center_point.z();

		if (debug_info.height_diff <= debug_info.z_threshold) {
			debug_info.reason = "below_tool_tip_threshold";
			return false;
		}

		if (debug_info.height_diff > tool.total_height) {
			debug_info.reason = "beyond_tool_total_height";
			return true;
		}

		if (!target_cell_vertices.empty()) {
			const double dx = target_cell_vertices[0](0, 0) - center_point.x();
			const double dy = target_cell_vertices[0](1, 0) - center_point.y();
			const double dist_xy_sq = dx * dx + dy * dy;

			if (debug_info.height_diff > tool.cylinder_height_threshold) {
				if (dist_xy_sq > tool.carriage_check_radius_sq) {
					debug_info.reason = "coarse_reject_carriage_xy";
					debug_info.dist_xy_sq = dist_xy_sq;
					debug_info.vertex_index = 0;
					return false;
				}
			}
			else {
				if (dist_xy_sq > tool.cylinder_check_radius_sq) {
					debug_info.reason = "coarse_reject_cylinder_xy";
					debug_info.dist_xy_sq = dist_xy_sq;
					debug_info.vertex_index = 0;
					return false;
				}
			}
		}

		for (int vertex_index = 0; vertex_index < static_cast<int>(target_cell_vertices.size()); ++vertex_index) {
			const auto& vertex = target_cell_vertices[vertex_index];
			const double diff_z = vertex(2, 0) - center_point.z();
			if (diff_z <= debug_info.z_threshold) {
				continue;
			}

			const double dx = vertex(0, 0) - center_point.x();
			const double dy = vertex(1, 0) - center_point.y();
			const double dist_xy_sq = dx * dx + dy * dy;

			if (diff_z <= tool.cylinder_height_threshold) {
				if (dist_xy_sq < tool.cylinder_r_sq) {
					debug_info.reason = "cylinder_radius_collision";
					debug_info.height_diff = diff_z;
					debug_info.dist_xy_sq = dist_xy_sq;
					debug_info.vertex_index = vertex_index;
					return true;
				}
			}
			else if (diff_z <= tool.total_height) {
				if (dist_xy_sq < tool.carriage_r_sq) {
					debug_info.reason = "carriage_radius_collision";
					debug_info.height_diff = diff_z;
					debug_info.dist_xy_sq = dist_xy_sq;
					debug_info.vertex_index = vertex_index;
					return true;
				}
			}
		}

		debug_info.reason = "no_vertex_collision";
		return false;
	}

	void WriteObjColoredVertex(
		std::ofstream& obj,
		const Eigen::Vector3d& point,
		const std::array<double, 3>& color,
		int& next_vertex_index)
	{
		obj << "v " << point.x() << " " << point.y() << " " << point.z() << " "
			<< color[0] << " " << color[1] << " " << color[2] << "\n";
		++next_vertex_index;
	}

	void WriteObjLine(std::ofstream& obj, int a, int b)
	{
		obj << "l " << a << " " << b << "\n";
	}

	std::vector<int> WriteObjRing(
		std::ofstream& obj,
		int& next_vertex_index,
		const Eigen::Matrix3d& rotated_to_world,
		const Eigen::Vector3d& center_rotated,
		double z_offset,
		double radius,
		const std::array<double, 3>& color,
		int segment_count = 32)
	{
		std::vector<int> ring_indices;
		ring_indices.reserve(segment_count);
		const double kPi = 3.14159265358979323846;
		for (int i = 0; i < segment_count; ++i) {
			const double theta = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(segment_count);
			const Eigen::Vector3d point_rotated(
				center_rotated.x() + radius * std::cos(theta),
				center_rotated.y() + radius * std::sin(theta),
				center_rotated.z() + z_offset);
			ring_indices.push_back(next_vertex_index);
			WriteObjColoredVertex(obj, rotated_to_world * point_rotated, color, next_vertex_index);
		}
		for (int i = 0; i < segment_count; ++i) {
			WriteObjLine(obj, ring_indices[i], ring_indices[(i + 1) % segment_count]);
		}
		return ring_indices;
	}

	void WriteToolCollisionCaseObj(
		std::ofstream& obj,
		int& next_vertex_index,
		const Eigen::Matrix3d& rotated_to_world,
		const Eigen::Vector3d& center_rotated,
		const std::vector<Eigen::Vector3d>& blocker_cell_vertices_rotated,
		const Eigen::Vector3d& target_point_world,
		const cutter& tool,
		int case_index,
		int s_id,
		int cell_id,
		int orientation_id,
		int blocker_cell_id,
		const ToolCollisionDebugInfo& collision_info)
	{
		if (!obj.is_open()) {
			return;
		}

		obj << "g collision_case_" << case_index
			<< "_s_" << s_id
			<< "_cell_" << cell_id
			<< "_ori_" << orientation_id
			<< "_blocker_" << blocker_cell_id
			<< "\n";
		obj << "# reason " << collision_info.reason
			<< " height_diff " << collision_info.height_diff
			<< " dist_xy_sq " << collision_info.dist_xy_sq << "\n";

		const std::array<double, 3> axis_color = { 0.0, 0.85, 1.0 };
		const std::array<double, 3> cylinder_color = { 1.0, 0.55, 0.05 };
		const std::array<double, 3> carriage_color = { 1.0, 0.9, 0.05 };
		const std::array<double, 3> blocker_color = { 1.0, 0.05, 0.05 };
		const std::array<double, 3> target_color = { 1.0, 0.0, 0.9 };

		const int axis_start = next_vertex_index;
		WriteObjColoredVertex(obj, rotated_to_world * center_rotated, axis_color, next_vertex_index);
		const int axis_end = next_vertex_index;
		WriteObjColoredVertex(
			obj,
			rotated_to_world * (center_rotated + Eigen::Vector3d(0.0, 0.0, tool.total_height)),
			axis_color,
			next_vertex_index);
		WriteObjLine(obj, axis_start, axis_end);

		const auto cylinder_bottom = WriteObjRing(
			obj, next_vertex_index, rotated_to_world, center_rotated, 0.0, tool.cylinder_r, cylinder_color);
		const auto cylinder_top = WriteObjRing(
			obj, next_vertex_index, rotated_to_world, center_rotated, tool.cylinder_height_threshold, tool.cylinder_r, cylinder_color);
		const auto carriage_bottom = WriteObjRing(
			obj, next_vertex_index, rotated_to_world, center_rotated, tool.cylinder_height_threshold, tool.carriage_r, carriage_color);
		const auto carriage_top = WriteObjRing(
			obj, next_vertex_index, rotated_to_world, center_rotated, tool.total_height, tool.carriage_r, carriage_color);
		for (int i = 0; i < static_cast<int>(cylinder_bottom.size()); i += 4) {
			WriteObjLine(obj, cylinder_bottom[i], cylinder_top[i]);
			WriteObjLine(obj, carriage_bottom[i], carriage_top[i]);
		}

		if (!blocker_cell_vertices_rotated.empty()) {
			std::vector<int> blocker_indices;
			blocker_indices.reserve(blocker_cell_vertices_rotated.size());
			for (const auto& point_rotated : blocker_cell_vertices_rotated) {
				blocker_indices.push_back(next_vertex_index);
				WriteObjColoredVertex(obj, rotated_to_world * point_rotated, blocker_color, next_vertex_index);
			}
			for (int i = 0; i < static_cast<int>(blocker_indices.size()); ++i) {
				WriteObjLine(obj, blocker_indices[i], blocker_indices[(i + 1) % blocker_indices.size()]);
			}
		}

		const double target_marker_radius = std::max(0.25, tool.cylinder_r * 0.15);
		const int target_x0 = next_vertex_index;
		WriteObjColoredVertex(obj, target_point_world - Eigen::Vector3d(target_marker_radius, 0.0, 0.0), target_color, next_vertex_index);
		const int target_x1 = next_vertex_index;
		WriteObjColoredVertex(obj, target_point_world + Eigen::Vector3d(target_marker_radius, 0.0, 0.0), target_color, next_vertex_index);
		const int target_y0 = next_vertex_index;
		WriteObjColoredVertex(obj, target_point_world - Eigen::Vector3d(0.0, target_marker_radius, 0.0), target_color, next_vertex_index);
		const int target_y1 = next_vertex_index;
		WriteObjColoredVertex(obj, target_point_world + Eigen::Vector3d(0.0, target_marker_radius, 0.0), target_color, next_vertex_index);
		const int target_z0 = next_vertex_index;
		WriteObjColoredVertex(obj, target_point_world - Eigen::Vector3d(0.0, 0.0, target_marker_radius), target_color, next_vertex_index);
		const int target_z1 = next_vertex_index;
		WriteObjColoredVertex(obj, target_point_world + Eigen::Vector3d(0.0, 0.0, target_marker_radius), target_color, next_vertex_index);
		WriteObjLine(obj, target_x0, target_x1);
		WriteObjLine(obj, target_y0, target_y1);
		WriteObjLine(obj, target_z0, target_z1);
	}

	bool FindHighestZMappedVertex(
		const Eigen::MatrixXd& original_vertices,
		const std::vector<vasco::VoronoiCell>& voronoi_cells,
		const std::unordered_map<int, int>& map_S_and_vertex,
		const std::vector<bool>* active_cell_mask,
		const std::vector<bool>* searched_s_flags,
		int& s_id,
		int& cell_id)
	{
		s_id = -1;
		cell_id = -1;
		double max_z = MIN_D;
		for (const auto& entry : map_S_and_vertex) {
			const int candidate_s_id = entry.first;
			const int candidate_cell_id = entry.second;
			if (searched_s_flags != nullptr
				&& candidate_s_id >= 0
				&& candidate_s_id < static_cast<int>(searched_s_flags->size())
				&& (*searched_s_flags)[candidate_s_id]) {
				continue;
			}
			if (candidate_cell_id < 0
				|| candidate_cell_id >= original_vertices.rows()
				|| candidate_cell_id >= static_cast<int>(voronoi_cells.size())
				|| !voronoi_cells[candidate_cell_id].is_available
				|| voronoi_cells[candidate_cell_id].all_points_in_polygon.empty()) {
				continue;
			}
			if (active_cell_mask != nullptr
				&& candidate_cell_id < static_cast<int>(active_cell_mask->size())
				&& !(*active_cell_mask)[candidate_cell_id]) {
				continue;
			}

			const double z = original_vertices(candidate_cell_id, 2);
			if (z > max_z) {
				max_z = z;
				s_id = candidate_s_id;
				cell_id = candidate_cell_id;
			}
		}

		return s_id >= 0 && cell_id >= 0;
	}

	std::string MakeToolCollisionObjBasePath(const std::string& output_file)
	{
		const std::string extension = ".obj";
		if (output_file.size() >= extension.size()
			&& output_file.substr(output_file.size() - extension.size()) == extension) {
			return output_file.substr(0, output_file.size() - extension.size());
		}
		return output_file;
	}

	void WriteHighestZInaccessiblePointAllOrientationToolCollisionObj(
		const std::string& output_file,
		const Eigen::MatrixXd& original_vertices,
		const std::vector<vasco::VoronoiCell>& voronoi_cells,
		const std::vector<Eigen::Vector3d>& orientation_samples,
		const std::unordered_map<int, int>& map_S_and_vertex,
		cutter tool,
		int& selected_s_id,
		int& selected_cell_id,
		const std::vector<bool>* active_cell_mask,
		const std::vector<bool>* searched_s_flags)
	{
		PrepareToolForCollision(tool);
		selected_s_id = -1;
		selected_cell_id = -1;
		if (!FindHighestZMappedVertex(
			original_vertices,
			voronoi_cells,
			map_S_and_vertex,
			active_cell_mask,
			searched_s_flags,
			selected_s_id,
			selected_cell_id)) {
			return;
		}

		const std::string output_base = MakeToolCollisionObjBasePath(output_file);
		const std::string index_file = output_base + "_index.csv";
		std::ofstream index_report(index_file);
		if (index_report.is_open()) {
			index_report << std::setprecision(17);
			index_report << "orientation_id,orientation_x,orientation_y,orientation_z,obj_file,blocker_cell,reason,height_diff,dist_xy_sq,center_x,center_y,center_z\n";
		}

		const Eigen::Vector3d target_point_world(
			original_vertices(selected_cell_id, 0),
			original_vertices(selected_cell_id, 1),
			original_vertices(selected_cell_id, 2));

		for (int ori = 0; ori < static_cast<int>(orientation_samples.size()); ++ori) {
			const Eigen::Vector3d orientation_vector = orientation_samples[ori].normalized();
			const Eigen::Matrix3d rot_matrix =
				Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(0, 0, 1), orientation_vector).toRotationMatrix();
			const Eigen::Matrix3d rot_matrix_inverse = rot_matrix.inverse();

			std::vector<Eigen::Vector3d> rotated_sites(original_vertices.rows());
			for (int site_id = 0; site_id < original_vertices.rows(); ++site_id) {
				rotated_sites[site_id] = rot_matrix_inverse * Eigen::Vector3d(
					original_vertices(site_id, 0),
					original_vertices(site_id, 1),
					original_vertices(site_id, 2));
			}

			std::vector<std::vector<Eigen::Vector3d>> rotated_cell_vertices(voronoi_cells.size());
			std::vector<double> max_z_of_cells(voronoi_cells.size(), MIN_D);
			for (int other_cell_id = 0; other_cell_id < static_cast<int>(voronoi_cells.size()); ++other_cell_id) {
				rotated_cell_vertices[other_cell_id].reserve(voronoi_cells[other_cell_id].all_points_in_polygon.size());
				for (const auto& point : voronoi_cells[other_cell_id].all_points_in_polygon) {
					Eigen::Vector3d rotated_point = rot_matrix_inverse * Eigen::Vector3d(point.x(), point.y(), point.z());
					rotated_cell_vertices[other_cell_id].push_back(rotated_point);
					max_z_of_cells[other_cell_id] = std::max(max_z_of_cells[other_cell_id], rotated_point.z());
				}
			}

			Eigen::Vector3d normal(0.0, 0.0, 0.0);
			const auto& selected_cell_vertices = rotated_cell_vertices[selected_cell_id];
			for (int j = 0; j < static_cast<int>(selected_cell_vertices.size()); ++j) {
				const Eigen::Vector3d& v1 = rotated_sites[voronoi_cells[selected_cell_id].site];
				const Eigen::Vector3d& v2 = selected_cell_vertices[j];
				const Eigen::Vector3d& v3 = selected_cell_vertices[(j + 1) % selected_cell_vertices.size()];
				normal += (v2 - v1).cross(v3 - v1);
			}
			if (normal.norm() <= 1e-12) {
				continue;
			}
			normal.normalize();

			const Eigen::Vector3d center_point = ComputeToolCenter(
				rotated_sites[voronoi_cells[selected_cell_id].site],
				normal,
				tool.cylinder_r);

			int blocker_cell = -1;
			ToolCollisionDebugInfo collision_info;
			for (int other_cell_id = 0; other_cell_id < static_cast<int>(voronoi_cells.size()); ++other_cell_id) {
				if (selected_cell_id == other_cell_id || !voronoi_cells[other_cell_id].is_available) {
					continue;
				}
				if (active_cell_mask != nullptr
					&& other_cell_id < static_cast<int>(active_cell_mask->size())
					&& !(*active_cell_mask)[other_cell_id]) {
					continue;
				}
				if (CheckToolCollisionWithCellForDebug(
					center_point,
					rotated_cell_vertices[other_cell_id],
					max_z_of_cells[other_cell_id],
					tool,
					collision_info)) {
					blocker_cell = other_cell_id;
					break;
				}
			}

			std::vector<Eigen::Vector3d> blocker_vertices;
			if (blocker_cell >= 0) {
				blocker_vertices = rotated_cell_vertices[blocker_cell];
			}
			else {
				collision_info.reason = "no_collision";
			}

			const std::string orientation_obj_file = output_base + "_ori_" + std::to_string(ori) + ".obj";
			std::ofstream obj(orientation_obj_file);
			if (!obj.is_open()) {
				std::cout << "[AccessibilityDebug] cannot open file for writing: " << orientation_obj_file << std::endl;
				continue;
			}

			obj << std::setprecision(17);
			obj << "# Tool collision visualization for the highest-z all-direction-inaccessible point\n";
			obj << "# selected_s_id " << selected_s_id << "\n";
			obj << "# selected_cell_id " << selected_cell_id << "\n";
			obj << "# selected_point "
				<< target_point_world.x() << " "
				<< target_point_world.y() << " "
				<< target_point_world.z() << "\n";
			obj << "# orientation_id " << ori << "\n";
			obj << "# orientation_vector "
				<< orientation_vector.x() << " "
				<< orientation_vector.y() << " "
				<< orientation_vector.z() << "\n";
			obj << "# Cyan line: tool axis\n";
			obj << "# Green line: sampled tool direction vector from the selected point\n";
			obj << "# Orange rings: cylinder radius envelope\n";
			obj << "# Yellow ring: carriage radius envelope\n";
			obj << "# Red loop: first blocker Voronoi cell for this orientation, if any\n";
			obj << "# Magenta cross: selected inaccessible point\n";

			int next_vertex_index = 1;
			const std::array<double, 3> direction_color = { 0.0, 1.0, 0.25 };
			const int direction_start = next_vertex_index;
			WriteObjColoredVertex(obj, target_point_world, direction_color, next_vertex_index);
			const int direction_end = next_vertex_index;
			WriteObjColoredVertex(
				obj,
				target_point_world + orientation_vector * std::max(tool.total_height, tool.cylinder_r * 3.0),
				direction_color,
				next_vertex_index);
			WriteObjLine(obj, direction_start, direction_end);

			WriteToolCollisionCaseObj(
				obj,
				next_vertex_index,
				rot_matrix,
				center_point,
				blocker_vertices,
				target_point_world,
				tool,
				ori,
				selected_s_id,
				selected_cell_id,
				ori,
				blocker_cell,
				collision_info);

			if (index_report.is_open()) {
				index_report << ori << ","
					<< orientation_vector.x() << ","
					<< orientation_vector.y() << ","
					<< orientation_vector.z() << ","
					<< orientation_obj_file << ","
					<< blocker_cell << ","
					<< collision_info.reason << ","
					<< collision_info.height_diff << ","
					<< collision_info.dist_xy_sq << ","
					<< center_point.x() << ","
					<< center_point.y() << ","
					<< center_point.z() << "\n";
			}
		}

		std::cout << "[AccessibilityDebug] wrote highest-z per-orientation tool collision visualizations with base "
			<< output_base << "_ori_<id>.obj" << std::endl;
	}

	void WriteSubtractiveAllDirectionReasonDiagnostics(
		const std::string& vis_dir,
		const std::string& file_token,
		const Eigen::MatrixXd& original_vertices,
		const std::vector<vasco::VoronoiCell>& voronoi_cells,
		const std::vector<Eigen::Vector3d>& orientation_samples,
		const std::unordered_map<int, int>& map_S_and_vertex,
		cutter tool)
	{
		if (map_S_and_vertex.empty() || orientation_samples.empty()) {
			return;
		}

		PrepareToolForCollision(tool);
		const double horizontal_z_threshold = 0.2;
		std::vector<int> horizontal_orientation_ids;
		for (int ori = 0; ori < static_cast<int>(orientation_samples.size()); ++ori) {
			if (std::abs(orientation_samples[ori].z()) <= horizontal_z_threshold) {
				horizontal_orientation_ids.push_back(ori);
			}
		}
		if (horizontal_orientation_ids.empty()) {
			return;
		}

		const std::string output_file =
			JoinPath(vis_dir, "access_debug_subtractive_all_direction_reason_" + file_token + ".csv");
		const std::string meta_file =
			JoinPath(vis_dir, "access_debug_subtractive_all_direction_reason_" + file_token + "_meta.txt");
		const std::string tool_collision_obj_file =
			JoinPath(vis_dir, "access_debug_subtractive_initial_full_model_tool_collision_" + file_token + ".obj");
		EnsureParentDirectory(output_file);
		std::ofstream report(output_file);
		if (!report.is_open()) {
			std::cout << "[AccessibilityDebug] cannot open file for writing: " << output_file << std::endl;
			return;
		}
		std::ofstream meta_report(meta_file);

		report << std::setprecision(17);
		const int max_diagnostic_points = 500;
		int highest_z_s_id = -1;
		int highest_z_cell_id = -1;
		WriteHighestZInaccessiblePointAllOrientationToolCollisionObj(
			tool_collision_obj_file,
			original_vertices,
			voronoi_cells,
			orientation_samples,
			map_S_and_vertex,
			tool,
			highest_z_s_id,
			highest_z_cell_id);
		if (meta_report.is_open()) {
			meta_report << std::setprecision(17);
			meta_report << "horizontal_z_threshold: " << horizontal_z_threshold << "\n";
			meta_report << "horizontal_orientation_count: " << horizontal_orientation_ids.size() << "\n";
			meta_report << "max_diagnostic_points: " << max_diagnostic_points << "\n";
			meta_report << "tool_collision_obj_scope: initial full model, highest-z all-direction-inaccessible point, one OBJ per sampled orientation\n";
			meta_report << "tool_collision_obj_selected_s_id: " << highest_z_s_id << "\n";
			meta_report << "tool_collision_obj_selected_cell_id: " << highest_z_cell_id << "\n";
			meta_report << "tool_collision_obj_base: " << MakeToolCollisionObjBasePath(tool_collision_obj_file) << "_ori_<id>.obj\n";
			meta_report << "tool_collision_obj_index: " << MakeToolCollisionObjBasePath(tool_collision_obj_file) << "_index.csv\n";
			meta_report << "tool.cylinder_r: " << tool.cylinder_r << "\n";
			meta_report << "tool.cylinder_height: " << tool.cylinder_height << "\n";
			meta_report << "tool.ball_r: " << tool.ball_r << "\n";
			meta_report << "tool.carriage_r: " << tool.carriage_r << "\n";
			meta_report << "tool.carriage_height: " << tool.carriage_height << "\n";
			meta_report << "tool.total_height: " << tool.total_height << "\n";
		}
		report << "s_id,vertex_id,x,y,z,horizontal_accessible_count,"
			<< "beyond_total_height_count,cylinder_collision_count,carriage_collision_count,other_collision_count,"
			<< "most_horizontal_ori,orientation_x,orientation_y,orientation_z,first_blocker_cell,first_reason,"
			<< "first_height_diff,first_dist_xy_sq,first_center_x,first_center_y,first_center_z\n";

		std::vector<std::pair<int, int>> sorted_entries(map_S_and_vertex.begin(), map_S_and_vertex.end());
		std::sort(sorted_entries.begin(), sorted_entries.end());
		int diagnostic_point_count = 0;
		for (const auto& entry : sorted_entries) {
			if (diagnostic_point_count >= max_diagnostic_points) {
				break;
			}
			++diagnostic_point_count;

			const int s_id = entry.first;
			const int cell_id = entry.second;
			if (cell_id < 0 || cell_id >= static_cast<int>(voronoi_cells.size())
				|| cell_id >= original_vertices.rows()
				|| !voronoi_cells[cell_id].is_available
				|| voronoi_cells[cell_id].all_points_in_polygon.empty()) {
				continue;
			}

			int horizontal_accessible_count = 0;
			int beyond_total_height_count = 0;
			int cylinder_collision_count = 0;
			int carriage_collision_count = 0;
			int other_collision_count = 0;
			int most_horizontal_ori = -1;
			double most_horizontal_abs_z = MAX_D;
			int first_blocker_cell = -1;
			ToolCollisionDebugInfo first_collision_info;
			Eigen::Vector3d first_center(0.0, 0.0, 0.0);

			for (int ori : horizontal_orientation_ids) {
				if (std::abs(orientation_samples[ori].z()) < most_horizontal_abs_z) {
					most_horizontal_abs_z = std::abs(orientation_samples[ori].z());
					most_horizontal_ori = ori;
				}

				const Eigen::Matrix3d rot_matrix =
					Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(0, 0, 1), orientation_samples[ori]).toRotationMatrix();
				const Eigen::Matrix3d rot_matrix_inverse = rot_matrix.inverse();

				std::vector<Eigen::Vector3d> rotated_sites(original_vertices.rows());
				std::vector<std::vector<Eigen::Vector3d>> rotated_cell_vertices(voronoi_cells.size());
				std::vector<double> max_z_of_cells(voronoi_cells.size(), MIN_D);
				for (int site_id = 0; site_id < original_vertices.rows(); ++site_id) {
					rotated_sites[site_id] = rot_matrix_inverse * Eigen::Vector3d(
						original_vertices(site_id, 0),
						original_vertices(site_id, 1),
						original_vertices(site_id, 2));
				}
				for (int other_cell_id = 0; other_cell_id < static_cast<int>(voronoi_cells.size()); ++other_cell_id) {
					rotated_cell_vertices[other_cell_id].reserve(voronoi_cells[other_cell_id].all_points_in_polygon.size());
					for (const auto& point : voronoi_cells[other_cell_id].all_points_in_polygon) {
						Eigen::Vector3d rotated_point = rot_matrix_inverse * Eigen::Vector3d(point.x(), point.y(), point.z());
						rotated_cell_vertices[other_cell_id].push_back(rotated_point);
						max_z_of_cells[other_cell_id] = std::max(max_z_of_cells[other_cell_id], rotated_point.z());
					}
				}

				Eigen::Vector3d normal(0.0, 0.0, 0.0);
				const auto& cell_vertices = rotated_cell_vertices[cell_id];
				for (int j = 0; j < static_cast<int>(cell_vertices.size()); ++j) {
					const Eigen::Vector3d& v1 = rotated_sites[voronoi_cells[cell_id].site];
					const Eigen::Vector3d& v2 = cell_vertices[j];
					const Eigen::Vector3d& v3 = cell_vertices[(j + 1) % cell_vertices.size()];
					normal += (v2 - v1).cross(v3 - v1);
				}
				if (normal.norm() <= 1e-12) {
					++other_collision_count;
					continue;
				}
				normal.normalize();

				const Eigen::Vector3d center_point = ComputeToolCenter(
					rotated_sites[voronoi_cells[cell_id].site],
					normal,
					tool.cylinder_r);

				bool collision = false;
				ToolCollisionDebugInfo collision_info;
				int blocker_cell = -1;
				for (int other_cell_id = 0; other_cell_id < static_cast<int>(voronoi_cells.size()); ++other_cell_id) {
					if (cell_id == other_cell_id || !voronoi_cells[other_cell_id].is_available) {
						continue;
					}

					if (CheckToolCollisionWithCellForDebug(
						center_point,
						rotated_cell_vertices[other_cell_id],
						max_z_of_cells[other_cell_id],
						tool,
						collision_info)) {
						collision = true;
						blocker_cell = other_cell_id;
						break;
					}
				}

				if (!collision) {
					++horizontal_accessible_count;
					continue;
				}

				if (first_blocker_cell == -1 && ori == most_horizontal_ori) {
					first_blocker_cell = blocker_cell;
					first_collision_info = collision_info;
					first_center = center_point;
				}

				if (collision_info.reason == "beyond_tool_total_height") {
					++beyond_total_height_count;
				}
				else if (collision_info.reason == "cylinder_radius_collision") {
					++cylinder_collision_count;
				}
				else if (collision_info.reason == "carriage_radius_collision") {
					++carriage_collision_count;
				}
				else {
					++other_collision_count;
				}
			}

			if (most_horizontal_ori != -1) {
				first_blocker_cell = -1;
				first_collision_info = ToolCollisionDebugInfo{};
				first_center = Eigen::Vector3d(0.0, 0.0, 0.0);
				const int ori = most_horizontal_ori;
				const Eigen::Matrix3d rot_matrix =
					Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(0, 0, 1), orientation_samples[ori]).toRotationMatrix();
				const Eigen::Matrix3d rot_matrix_inverse = rot_matrix.inverse();
				std::vector<Eigen::Vector3d> rotated_sites(original_vertices.rows());
				std::vector<std::vector<Eigen::Vector3d>> rotated_cell_vertices(voronoi_cells.size());
				std::vector<double> max_z_of_cells(voronoi_cells.size(), MIN_D);
				for (int site_id = 0; site_id < original_vertices.rows(); ++site_id) {
					rotated_sites[site_id] = rot_matrix_inverse * Eigen::Vector3d(
						original_vertices(site_id, 0),
						original_vertices(site_id, 1),
						original_vertices(site_id, 2));
				}
				for (int other_cell_id = 0; other_cell_id < static_cast<int>(voronoi_cells.size()); ++other_cell_id) {
					for (const auto& point : voronoi_cells[other_cell_id].all_points_in_polygon) {
						Eigen::Vector3d rotated_point = rot_matrix_inverse * Eigen::Vector3d(point.x(), point.y(), point.z());
						rotated_cell_vertices[other_cell_id].push_back(rotated_point);
						max_z_of_cells[other_cell_id] = std::max(max_z_of_cells[other_cell_id], rotated_point.z());
					}
				}
				Eigen::Vector3d normal(0.0, 0.0, 0.0);
				const auto& cell_vertices = rotated_cell_vertices[cell_id];
				for (int j = 0; j < static_cast<int>(cell_vertices.size()); ++j) {
					const Eigen::Vector3d& v1 = rotated_sites[voronoi_cells[cell_id].site];
					const Eigen::Vector3d& v2 = cell_vertices[j];
					const Eigen::Vector3d& v3 = cell_vertices[(j + 1) % cell_vertices.size()];
					normal += (v2 - v1).cross(v3 - v1);
				}
				if (normal.norm() > 1e-12) {
					normal.normalize();
					first_center = ComputeToolCenter(rotated_sites[voronoi_cells[cell_id].site], normal, tool.cylinder_r);
					for (int other_cell_id = 0; other_cell_id < static_cast<int>(voronoi_cells.size()); ++other_cell_id) {
						if (cell_id == other_cell_id || !voronoi_cells[other_cell_id].is_available) {
							continue;
						}
						if (CheckToolCollisionWithCellForDebug(
							first_center,
							rotated_cell_vertices[other_cell_id],
							max_z_of_cells[other_cell_id],
							tool,
							first_collision_info)) {
							first_blocker_cell = other_cell_id;
							break;
						}
					}
				}
			}

			const Eigen::Vector3d ori_vec =
				(most_horizontal_ori >= 0) ? orientation_samples[most_horizontal_ori] : Eigen::Vector3d(0, 0, 0);
			report << s_id << "," << cell_id << ","
				<< original_vertices(cell_id, 0) << ","
				<< original_vertices(cell_id, 1) << ","
				<< original_vertices(cell_id, 2) << ","
				<< horizontal_accessible_count << ","
				<< beyond_total_height_count << ","
				<< cylinder_collision_count << ","
				<< carriage_collision_count << ","
				<< other_collision_count << ","
				<< most_horizontal_ori << ","
				<< ori_vec.x() << "," << ori_vec.y() << "," << ori_vec.z() << ","
				<< first_blocker_cell << ","
				<< first_collision_info.reason << ","
				<< first_collision_info.height_diff << ","
				<< first_collision_info.dist_xy_sq << ","
				<< first_center.x() << ","
				<< first_center.y() << ","
				<< first_center.z() << "\n";
		}

		std::cout << "[AccessibilityDebug] wrote subtractive all-direction reason report to "
			<< output_file << std::endl;
	}

	void WriteSubtractiveAccessibilityDebugVisualizations(
		const std::string& vis_dir,
		const std::string& node_tag,
		double marker_radius,
		const Eigen::MatrixXd& original_vertices,
		const std::vector<bool>& judge_S_be_searched,
		const std::vector<bool>& judge_covering_points_be_searched,
		const std::unordered_map<int, int>& map_S_and_vertex,
		const std::unordered_map<int, int>& map_covering_points_and_vertex)
	{
		std::vector<Eigen::Vector3d> remaining_subtractive_points;
		for (int point_id = 0; point_id < static_cast<int>(judge_S_be_searched.size()); ++point_id) {
			if (judge_S_be_searched[point_id]) {
				continue;
			}

			const auto map_it = map_S_and_vertex.find(point_id);
			if (map_it == map_S_and_vertex.end()) {
				continue;
			}

			const int vertex_id = map_it->second;
			if (vertex_id < 0 || vertex_id >= original_vertices.rows()) {
				continue;
			}

			remaining_subtractive_points.emplace_back(
				original_vertices(vertex_id, 0),
				original_vertices(vertex_id, 1),
				original_vertices(vertex_id, 2));
		}
		WriteDebugMarkersObj(
			JoinPath(vis_dir, "access_debug_subtractive_remaining_S" + node_tag + ".obj"),
			remaining_subtractive_points,
			{ 1.0, 0.05, 0.05 },
			marker_radius);
		WriteDebugMarkersObj(
			JoinPath(vis_dir, "access_debug_subtractive_all_direction_remaining" + node_tag + ".obj"),
			remaining_subtractive_points,
			{ 0.95, 0.05, 0.85 },
			marker_radius);

		std::vector<Eigen::Vector3d> remaining_covering_points;
		for (int point_id = 0; point_id < static_cast<int>(judge_covering_points_be_searched.size()); ++point_id) {
			if (judge_covering_points_be_searched[point_id]) {
				continue;
			}

			const auto map_it = map_covering_points_and_vertex.find(point_id);
			if (map_it == map_covering_points_and_vertex.end()) {
				continue;
			}

			const int vertex_id = map_it->second;
			if (vertex_id < 0 || vertex_id >= original_vertices.rows()) {
				continue;
			}

			remaining_covering_points.emplace_back(
				original_vertices(vertex_id, 0),
				original_vertices(vertex_id, 1),
				original_vertices(vertex_id, 2));
		}
		WriteDebugMarkersObj(
			JoinPath(vis_dir, "access_debug_subtractive_remaining_covering" + node_tag + ".obj"),
			remaining_covering_points,
			{ 0.05, 0.35, 1.0 },
			marker_radius);
	}

	std::vector<bool> BuildActiveOriginalVertexMask(
		const Eigen::MatrixXd& original_vertices,
		const std::vector<Vertex>& current_vertices,
		double eps)
	{
		std::vector<bool> active(original_vertices.rows(), false);
		for (int i = 0; i < original_vertices.rows(); ++i) {
			for (const auto& vertex : current_vertices) {
				if (std::abs(original_vertices(i, 0) - vertex.x) <= eps
					&& std::abs(original_vertices(i, 1) - vertex.y) <= eps
					&& std::abs(original_vertices(i, 2) - vertex.z) <= eps) {
					active[i] = true;
					break;
				}
			}
		}
		return active;
	}

	void WriteAdditiveAccessibilityDebugVisualization(
		const std::string& vis_dir,
		const Layer_Graph& layer_graph,
		const std::string& node_tag,
		int orientation_id,
		double marker_radius)
	{
		std::set<int> additive_collision_target_nodes;
		for (int edge_id = layer_graph.cont_normal_dependency_edges;
			edge_id < static_cast<int>(layer_graph.temp_edges.size());
			++edge_id) {
			additive_collision_target_nodes.insert(layer_graph.temp_edges[edge_id].second);
		}

		std::vector<Eigen::Vector3d> additive_collision_points;
		additive_collision_points.reserve(additive_collision_target_nodes.size());
		for (int node_id : additive_collision_target_nodes) {
			Eigen::Vector3d centroid;
			if (GetLayerGraphNodeCentroid(layer_graph, node_id, centroid)) {
				additive_collision_points.push_back(centroid);
			}
		}

		if (additive_collision_points.empty()) {
			return;
		}

		WriteDebugMarkersObj(
			JoinPath(vis_dir, "access_debug_additive_collision_targets"
				+ node_tag
				+ "_ori" + std::to_string(orientation_id) + ".obj"),
			additive_collision_points,
			{ 1.0, 0.55, 0.05 },
			marker_radius);
	}

	void WriteAdditiveRootLayerSelfSupportDebugVisualization(
		const std::string& vis_dir,
		const Layer_Graph& layer_graph,
		const std::string& node_tag,
		int orientation_id,
		const Eigen::Vector3d& orientation,
		const SurfaceMesh& mesh)
	{
		const int node_count = layer_graph.total_node_num;
		if (node_count <= 0) {
			return;
		}
		const Eigen::Matrix3d slicing_to_original =
			Eigen::Quaterniond::FromTwoVectors(
				Eigen::Vector3d(0.0, 0.0, 1.0),
				orientation).toRotationMatrix();

		const std::string file_stem =
			"access_debug_additive_root_self_support"
			+ node_tag
			+ "_ori" + std::to_string(orientation_id);
		const std::string obj_file = JoinPath(vis_dir, file_stem + ".obj");
		const std::string csv_file = JoinPath(vis_dir, file_stem + ".csv");
		const std::string downward_face_csv_file =
			JoinPath(vis_dir, file_stem + "_downward_faces.csv");
		EnsureParentDirectory(obj_file);

		double min_layer_z = std::numeric_limits<double>::infinity();
		double max_layer_z = -std::numeric_limits<double>::infinity();
		for (int node_id = 0; node_id < node_count; ++node_id) {
			const auto index_it = layer_graph.data.index.find(node_id);
			if (index_it == layer_graph.data.index.end()) {
				continue;
			}
			const int slice_id = index_it->second.first;
			const int component_id = index_it->second.second;
			if (slice_id < 0
				|| slice_id >= static_cast<int>(layer_graph.data.z_value.size())
				|| component_id < 0
				|| component_id >= static_cast<int>(layer_graph.data.z_value[slice_id].size())
				|| layer_graph.data.z_value[slice_id][component_id].empty()) {
				continue;
			}
			const double z = layer_graph.data.z_value[slice_id][component_id].front();
			min_layer_z = std::min(min_layer_z, z);
			max_layer_z = std::max(max_layer_z, z);
		}
		if (!std::isfinite(min_layer_z) || !std::isfinite(max_layer_z)) {
			return;
		}

		std::ofstream obj(obj_file);
		std::ofstream csv(csv_file);
		if (!obj.is_open() || !csv.is_open()) {
			std::cout << "[AccessibilityDebug] cannot open additive root-layer debug output: "
				<< file_stem << std::endl;
			return;
		}

		obj << std::setprecision(17);
		obj << "# Additive DFS root-layer self-support diagnostics\n";
		obj << "# OBJ coordinates have been transformed back to the original model frame.\n";
		obj << "# Layer z values and self-support normal tests remain in the slicing frame.\n";
		obj << "# Gray: minimum/maximum-z reference layer; green: self-supporting root candidate;\n";
		obj << "# red: production non-self-supporting root candidate allowed by the root exemption;\n";
		obj << "# magenta: production says self-supporting, but the counted normals exceed the downward-face threshold.\n";
		obj << "# Orange faces: too-downward faces associated only with non-root layers.\n";
		obj << "# Magenta faces: too-downward faces associated with at least one root candidate.\n";
		obj << "# orientation " << orientation.x() << " " << orientation.y() << " " << orientation.z() << "\n";
		obj << "# too_downward_normal_z_threshold "
			<< -std::sin(std::acos(-1.0) / 3.6) << "\n";
		obj << "# too_downward_face_ratio_threshold 0.2\n";
		obj << "# min_layer_z " << min_layer_z << "\n";
		obj << "# max_layer_z " << max_layer_z << "\n";

		csv << std::setprecision(17);
		csv << "node_id,slice_id,component_id,z,is_root_candidate,is_self_support,"
			"would_be_self_support_from_face_counts,orientation_x,orientation_y,orientation_z,"
			"initial_in_degree,initial_out_degree,point_count,hole_count,"
			"associated_face_count,stored_face_normal_count,"
			"valid_normal_count,invalid_face_count,"
			"too_downward_face_count,too_downward_face_ratio,"
			"too_downward_normal_z_threshold,too_downward_face_ratio_threshold,"
			"is_min_z,is_max_z\n";

		int next_vertex_id = 1;
		int root_candidate_count = 0;
		int exempt_candidate_count = 0;
		int debug_non_self_support_candidate_count = 0;
		const double z_eps = 1e-9 * std::max(
			1.0,
			std::max(std::abs(min_layer_z), std::abs(max_layer_z)));
		for (int node_id = 0; node_id < node_count; ++node_id) {
			const auto index_it = layer_graph.data.index.find(node_id);
			if (index_it == layer_graph.data.index.end()) {
				continue;
			}

			const int slice_id = index_it->second.first;
			const int component_id = index_it->second.second;
			if (slice_id < 0
				|| slice_id >= static_cast<int>(layer_graph.data.slice_points.size())
				|| slice_id >= static_cast<int>(layer_graph.data.slice_points_holes.size())
				|| slice_id >= static_cast<int>(layer_graph.data.z_value.size())
				|| component_id < 0
				|| component_id >= static_cast<int>(layer_graph.data.slice_points[slice_id].size())
				|| component_id >= static_cast<int>(layer_graph.data.slice_points_holes[slice_id].size())
				|| component_id >= static_cast<int>(layer_graph.data.z_value[slice_id].size())
				|| layer_graph.data.z_value[slice_id][component_id].empty()) {
				continue;
			}

			const bool is_root_candidate =
				node_id < static_cast<int>(layer_graph.out_degree.size())
				&& layer_graph.out_degree[node_id] == 0;
			const bool is_self_support =
				node_id < static_cast<int>(layer_graph.is_the_layer_self_suppot.size())
				&& layer_graph.is_the_layer_self_suppot[node_id];
			const double z = layer_graph.data.z_value[slice_id][component_id].front();
			const bool is_min_z = std::abs(z - min_layer_z) <= z_eps;
			const bool is_max_z = std::abs(z - max_layer_z) <= z_eps;
			const auto& outer = layer_graph.data.slice_points[slice_id][component_id];
			const auto& holes = layer_graph.data.slice_points_holes[slice_id][component_id];
			const int in_degree =
				node_id < static_cast<int>(layer_graph.in_degree.size())
				? layer_graph.in_degree[node_id]
				: -1;
			const int out_degree =
				node_id < static_cast<int>(layer_graph.out_degree.size())
				? layer_graph.out_degree[node_id]
				: -1;
			const int associated_face_count =
				node_id < static_cast<int>(layer_graph.self_support_associated_face_count.size())
				? layer_graph.self_support_associated_face_count[node_id]
				: 0;
			const int stored_normal_count =
				node_id < static_cast<int>(layer_graph.self_support_stored_normal_count.size())
				? layer_graph.self_support_stored_normal_count[node_id]
				: 0;
			const int valid_normal_count =
				node_id < static_cast<int>(layer_graph.self_support_valid_normal_count.size())
				? layer_graph.self_support_valid_normal_count[node_id]
				: 0;
			const int invalid_face_count =
				node_id < static_cast<int>(layer_graph.self_support_invalid_face_count.size())
				? layer_graph.self_support_invalid_face_count[node_id]
				: 0;
			const int too_downward_face_count =
				node_id < static_cast<int>(layer_graph.self_support_too_downward_face_count.size())
				? layer_graph.self_support_too_downward_face_count[node_id]
				: 0;
			const double too_downward_face_ratio =
				associated_face_count > 0
				? static_cast<double>(too_downward_face_count)
					/ static_cast<double>(associated_face_count)
				: 0.0;
			const double too_downward_normal_z_threshold =
				-std::sin(std::acos(-1.0) / 3.6);
			const bool would_be_self_support_from_face_counts =
				associated_face_count == 0 || too_downward_face_ratio <= 0.2;

			std::array<double, 3> color{ 0.55, 0.55, 0.55 };
			if (is_root_candidate) {
				++root_candidate_count;
				if (!is_self_support) {
					color = { 1.00, 0.08, 0.08 };
					++exempt_candidate_count;
				}
				else if (!would_be_self_support_from_face_counts) {
					color = { 1.00, 0.05, 0.85 };
					++debug_non_self_support_candidate_count;
				}
				else {
					color = { 0.10, 0.85, 0.20 };
				}
			}

			csv << node_id << ","
				<< slice_id << ","
				<< component_id << ","
				<< z << ","
				<< (is_root_candidate ? 1 : 0) << ","
				<< (is_self_support ? 1 : 0) << ","
				<< (would_be_self_support_from_face_counts ? 1 : 0) << ","
				<< orientation.x() << ","
				<< orientation.y() << ","
				<< orientation.z() << ","
				<< in_degree << ","
				<< out_degree << ","
				<< outer.size() << ","
				<< holes.size() << ","
				<< associated_face_count << ","
				<< stored_normal_count << ","
				<< valid_normal_count << ","
				<< invalid_face_count << ","
				<< too_downward_face_count << ","
				<< too_downward_face_ratio << ","
				<< too_downward_normal_z_threshold << ","
				<< 0.2 << ","
				<< (is_min_z ? 1 : 0) << ","
				<< (is_max_z ? 1 : 0) << "\n";

			if (!is_root_candidate && !is_min_z && !is_max_z) {
				continue;
			}

			auto write_loop = [&](const std::vector<Eigen::Vector2d>& loop, const std::string& group_name) {
				if (loop.size() < 2) {
					return;
				}

				obj << "g " << group_name << "\n";
				const int first_vertex_id = next_vertex_id;
				for (const auto& point : loop) {
					const Eigen::Vector3d original_point =
						slicing_to_original * Eigen::Vector3d(point.x(), point.y(), z);
					obj << "v "
						<< original_point.x() << " "
						<< original_point.y() << " "
						<< original_point.z() << " "
						<< color[0] << " " << color[1] << " " << color[2] << "\n";
					++next_vertex_id;
				}
				obj << "l";
				for (int vertex_id = first_vertex_id; vertex_id < next_vertex_id; ++vertex_id) {
					obj << " " << vertex_id;
				}
				obj << " " << first_vertex_id << "\n";
			};

			const std::string node_group =
				"node_" + std::to_string(node_id)
				+ "_slice_" + std::to_string(slice_id)
				+ (is_root_candidate
					? (!is_self_support
						? "_root_exempt"
						: (!would_be_self_support_from_face_counts
							? "_root_debug_downward_mismatch"
							: "_root_supported"))
					: (is_min_z ? "_min_z_reference" : "_max_z_reference"));
			write_loop(outer, node_group + "_outer");
			for (std::size_t hole_id = 0; hole_id < holes.size(); ++hole_id) {
				write_loop(holes[hole_id], node_group + "_hole_" + std::to_string(hole_id));
			}
		}

		// Collect each downward face once for geometry output while retaining
		// whether any root candidate references it for color classification.
		std::map<int, bool> downward_face_has_root_candidate;
		for (int node_id = 0; node_id < node_count; ++node_id) {
			if (node_id >= static_cast<int>(
				layer_graph.self_support_too_downward_face_ids.size())) {
				continue;
			}
			const bool is_root_candidate =
				node_id < static_cast<int>(layer_graph.out_degree.size())
				&& layer_graph.out_degree[node_id] == 0;
			for (int face_id : layer_graph.self_support_too_downward_face_ids[node_id]) {
				auto insertion = downward_face_has_root_candidate.emplace(
					face_id,
					is_root_candidate);
				if (!insertion.second && is_root_candidate) {
					insertion.first->second = true;
				}
			}
		}

		auto get_face_geometry = [&](int face_id,
			std::array<Eigen::Vector3d, 3>& points,
			Eigen::Vector3d& normal) {
			if (face_id < 0
				|| static_cast<std::size_t>(face_id) >= mesh.number_of_faces()) {
				return false;
			}
			const SurfaceMesh::Face_index face_index(face_id);
			const auto h0 = mesh.halfedge(face_index);
			const auto h1 = mesh.next(h0);
			const auto h2 = mesh.next(h1);
			const std::array<SurfaceMesh::Halfedge_index, 3> halfedges{ h0, h1, h2 };
			for (int i = 0; i < 3; ++i) {
				const Point_3& point = mesh.point(mesh.target(halfedges[i]));
				points[i] = Eigen::Vector3d(
					CGAL::to_double(point.x()),
					CGAL::to_double(point.y()),
					CGAL::to_double(point.z()));
			}
			normal = (points[1] - points[0]).cross(points[2] - points[0]);
			if (!normal.allFinite() || normal.squaredNorm() <= 1e-24) {
				return false;
			}
			normal.normalize();
			return true;
		};

		std::ofstream downward_face_csv(downward_face_csv_file);
		if (downward_face_csv.is_open()) {
			downward_face_csv << std::setprecision(17);
			downward_face_csv
				<< "node_id,slice_id,component_id,layer_z,is_root_candidate,"
				"face_id,slicing_normal_x,slicing_normal_y,slicing_normal_z,"
				"original_normal_x,original_normal_y,original_normal_z\n";
			for (int node_id = 0; node_id < node_count; ++node_id) {
				if (node_id >= static_cast<int>(
					layer_graph.self_support_too_downward_face_ids.size())) {
					continue;
				}
				const auto index_it = layer_graph.data.index.find(node_id);
				if (index_it == layer_graph.data.index.end()) {
					continue;
				}
				const int slice_id = index_it->second.first;
				const int component_id = index_it->second.second;
				if (slice_id < 0
					|| slice_id >= static_cast<int>(layer_graph.data.z_value.size())
					|| component_id < 0
					|| component_id >= static_cast<int>(
						layer_graph.data.z_value[slice_id].size())
					|| layer_graph.data.z_value[slice_id][component_id].empty()) {
					continue;
				}
				const bool is_root_candidate =
					node_id < static_cast<int>(layer_graph.out_degree.size())
					&& layer_graph.out_degree[node_id] == 0;
				const double layer_z =
					layer_graph.data.z_value[slice_id][component_id].front();
				for (int face_id :
					layer_graph.self_support_too_downward_face_ids[node_id]) {
					std::array<Eigen::Vector3d, 3> points;
					Eigen::Vector3d normal;
					if (!get_face_geometry(face_id, points, normal)) {
						continue;
					}
					const Eigen::Vector3d original_normal =
						slicing_to_original * normal;
					downward_face_csv
						<< node_id << ","
						<< slice_id << ","
						<< component_id << ","
						<< layer_z << ","
						<< (is_root_candidate ? 1 : 0) << ","
						<< face_id << ","
						<< normal.x() << ","
						<< normal.y() << ","
						<< normal.z() << ","
						<< original_normal.x() << ","
						<< original_normal.y() << ","
						<< original_normal.z() << "\n";
				}
			}
		}

		int written_downward_face_count = 0;
		for (const auto& face_entry : downward_face_has_root_candidate) {
			std::array<Eigen::Vector3d, 3> points;
			Eigen::Vector3d normal;
			if (!get_face_geometry(face_entry.first, points, normal)) {
				continue;
			}

			const std::array<double, 3> color = face_entry.second
				? std::array<double, 3>{ 1.00, 0.05, 0.85 }
				: std::array<double, 3>{ 1.00, 0.55, 0.05 };
			obj << "g too_downward_face_" << face_entry.first
				<< (face_entry.second ? "_root" : "_non_root") << "\n";
			for (const auto& point : points) {
				const Eigen::Vector3d original_point =
					slicing_to_original * point;
				obj << "v "
					<< original_point.x() << " "
					<< original_point.y() << " "
					<< original_point.z() << " "
					<< color[0] << " " << color[1] << " " << color[2] << "\n";
			}
			obj << "f " << next_vertex_id
				<< " " << next_vertex_id + 1
				<< " " << next_vertex_id + 2 << "\n";
			next_vertex_id += 3;
			++written_downward_face_count;
		}

		std::cout << "[AccessibilityDebug] additive DFS root-layer self-support diagnostics saved: "
			<< obj_file
			<< ", root_candidates=" << root_candidate_count
			<< ", exempt_candidates=" << exempt_candidate_count
			<< ", debug_non_self_support_candidates=" << debug_non_self_support_candidate_count
			<< ", too_downward_faces=" << written_downward_face_count
			<< std::endl;
	}

	void ShowMergedPatchGraphCutLabels(
		const std::string& vis_dir,
		const std::string& output_prefix,
		int height_of_beam_search,
		int orientation_count,
		const vasco::Slicer& merged_patch,
		const std::vector<int>& result_labels,
		const std::vector<Eigen::Vector3d>& orientation_samples)
	{
		if (orientation_count <= 0
			|| merged_patch.positions.empty()
			|| merged_patch.triangles.empty()) {
			std::cout << "[AccessibilityVisualization] Skip graph-cut Polyscope visualization: "
				<< "empty mesh or orientation set." << std::endl;
			return;
		}
		if (result_labels.size() != merged_patch.triangles.size()) {
			std::cerr << "[AccessibilityVisualization] Skip graph-cut Polyscope visualization: "
				<< "label_count=" << result_labels.size()
				<< " differs from face_count=" << merged_patch.triangles.size()
				<< "." << std::endl;
			return;
		}

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
		auto color_for_orientation = [](int orientation_id) {
			const std::size_t color_id = static_cast<std::size_t>(
				std::max(0, orientation_id)) % palette.size();
			return palette[color_id];
		};

		std::map<int, std::vector<vasco::Slicer::Tri3>> label_faces;
		int invalid_face_count = 0;
		for (std::size_t face_id = 0; face_id < merged_patch.triangles.size(); ++face_id) {
			const auto& triangle = merged_patch.triangles[face_id];
			bool valid = true;
			for (int vertex_id : triangle) {
				if (vertex_id < 0
					|| vertex_id >= static_cast<int>(merged_patch.positions.size())) {
					valid = false;
					break;
				}
				const auto& point = merged_patch.positions[vertex_id];
				if (!std::isfinite(point[0])
					|| !std::isfinite(point[1])
					|| !std::isfinite(point[2])) {
					valid = false;
					break;
				}
			}
			if (!valid) {
				++invalid_face_count;
				continue;
			}
			label_faces[result_labels[face_id]].push_back(triangle);
		}

		Eigen::Vector3d bounds_min(
			std::numeric_limits<double>::max(),
			std::numeric_limits<double>::max(),
			std::numeric_limits<double>::max());
		Eigen::Vector3d bounds_max(
			std::numeric_limits<double>::lowest(),
			std::numeric_limits<double>::lowest(),
			std::numeric_limits<double>::lowest());
		for (const auto& point : merged_patch.positions) {
			bounds_min.x() = std::min(bounds_min.x(), point[0]);
			bounds_min.y() = std::min(bounds_min.y(), point[1]);
			bounds_min.z() = std::min(bounds_min.z(), point[2]);
			bounds_max.x() = std::max(bounds_max.x(), point[0]);
			bounds_max.y() = std::max(bounds_max.y(), point[1]);
			bounds_max.z() = std::max(bounds_max.z(), point[2]);
		}
		const double model_diagonal = (bounds_max - bounds_min).norm();
		const double arrow_start_offset = std::max(1.0, 0.06 * model_diagonal);
		const double arrow_length = std::max(2.0, 0.12 * model_diagonal);

		std::vector<Eigen::Vector3d> arrow_points;
		std::vector<Eigen::Vector3d> arrow_directions;
		std::vector<std::vector<double>> arrow_colors;

		polyscope::removeAllStructures();
		for (const auto& label_entry : label_faces) {
			const int label_id = label_entry.first;
			const int patch_id = label_id / orientation_count + 1;
			const int orientation_id = label_id % orientation_count;
			const auto& faces = label_entry.second;
			if (faces.empty()) {
				continue;
			}

			const auto color = color_for_orientation(orientation_id);
			const std::string name_suffix =
				std::to_string(height_of_beam_search) + "_"
				+ std::to_string(patch_id) + "_"
				+ std::to_string(orientation_id);
			auto* mesh = polyscope::registerSurfaceMesh(
				output_prefix + "merged_patch_graphcut_patch_" + name_suffix,
				merged_patch.positions,
				faces);
			mesh->setSurfaceColor({
				static_cast<float>(color[0]),
				static_cast<float>(color[1]),
				static_cast<float>(color[2])
				});

			Eigen::Vector3d center = Eigen::Vector3d::Zero();
			std::size_t center_sample_count = 0;
			for (const auto& triangle : faces) {
				for (int vertex_id : triangle) {
					const auto& point = merged_patch.positions[vertex_id];
					center += Eigen::Vector3d(point[0], point[1], point[2]);
					++center_sample_count;
				}
			}
			if (center_sample_count == 0) {
				continue;
			}
			center /= static_cast<double>(center_sample_count);

			Eigen::Vector3d direction(0.0, 0.0, 1.0);
			if (orientation_id >= 0
				&& orientation_id < static_cast<int>(orientation_samples.size())
				&& orientation_samples[orientation_id].allFinite()
				&& orientation_samples[orientation_id].squaredNorm() > 1e-24) {
				direction = orientation_samples[orientation_id].normalized();
			}

			const Eigen::Vector3d arrow_start =
				center + direction * arrow_start_offset;
			const Eigen::Vector3d arrow_end =
				arrow_start + direction * arrow_length;
			Eigen::Vector3d reference_axis = std::abs(direction.z()) < 0.9
				? Eigen::Vector3d(0.0, 0.0, 1.0)
				: Eigen::Vector3d(0.0, 1.0, 0.0);
			Eigen::Vector3d perpendicular = direction.cross(reference_axis);
			if (perpendicular.squaredNorm() <= 1e-24) {
				perpendicular = Eigen::Vector3d(1.0, 0.0, 0.0);
			}
			else {
				perpendicular.normalize();
			}

			const double shaft_width = std::max(0.2, arrow_length * 0.08);
			const double head_length = arrow_length * 0.35;
			const double head_width = shaft_width * 2.0;
			const Eigen::Vector3d head_base =
				arrow_end - direction * head_length;
			const Eigen::Vector3d shaft_left =
				arrow_start + perpendicular * shaft_width;
			const Eigen::Vector3d shaft_right =
				arrow_start - perpendicular * shaft_width;
			const Eigen::Vector3d shaft_end_left =
				head_base + perpendicular * shaft_width;
			const Eigen::Vector3d shaft_end_right =
				head_base - perpendicular * shaft_width;
			const Eigen::Vector3d head_left =
				head_base + perpendicular * head_width;
			const Eigen::Vector3d head_right =
				head_base - perpendicular * head_width;

			std::vector<vasco::Slicer::Vec3> arrow_positions = {
				{shaft_left.x(), shaft_left.y(), shaft_left.z()},
				{shaft_right.x(), shaft_right.y(), shaft_right.z()},
				{shaft_end_right.x(), shaft_end_right.y(), shaft_end_right.z()},
				{shaft_end_left.x(), shaft_end_left.y(), shaft_end_left.z()},
				{head_left.x(), head_left.y(), head_left.z()},
				{head_right.x(), head_right.y(), head_right.z()},
				{arrow_end.x(), arrow_end.y(), arrow_end.z()}
			};
			std::vector<vasco::Slicer::Tri3> arrow_triangles = {
				{0, 1, 2}, {0, 2, 3},
				{4, 5, 6}, {3, 4, 6}, {2, 6, 5}
			};
			auto* arrow = polyscope::registerSurfaceMesh(
				output_prefix + "merged_patch_graphcut_arrow_" + name_suffix,
				arrow_positions,
				arrow_triangles);
			arrow->setSurfaceColor({
				static_cast<float>(color[0]),
				static_cast<float>(color[1]),
				static_cast<float>(color[2])
				});

			arrow_points.push_back(arrow_start);
			arrow_directions.push_back(direction);
			arrow_colors.push_back({color[0], color[1], color[2]});
		}

		if (!arrow_points.empty()) {
			Visual visual;
			visual.generateModelForRendering_11(
				arrow_points,
				arrow_directions,
				arrow_colors,
				JoinPath(
					vis_dir,
					output_prefix + "merged_patch_graphcut_patch_ori_arrows.obj"));
		}

		std::cout << "[AccessibilityVisualization] Showing graph-cut result in Polyscope: "
			<< "labels=" << label_faces.size()
			<< ", skipped_invalid_faces=" << invalid_face_count
			<< "." << std::endl;
		polyscope::show();
	}

}
