#include "data.h"
#include "vasco/core/ContourContainment.h"

void Data::ReadData(
	vector<vector<vector<Eigen::Vector3d>>> all_slice_points,
	vector<vector<vector<Eigen::Vector3d>>> all_slice_points_contain,
	const vector<vector<int>>& all_contour_parent_ids)
{
	const int total_layer = static_cast<int>(all_slice_points.size());
	slice_points.assign(total_layer, {});
	slice_points_contain.assign(total_layer, {});
	slice_points_holes.assign(total_layer, {});
	source_contour_ids.assign(total_layer, {});
	source_hole_contour_ids.assign(total_layer, {});
	contour_parent_ids.assign(total_layer, {});
	z_value.assign(total_layer, {});
	is_contour.assign(total_layer, {});
	adjacent_points.assign(total_layer, {});
	index.clear();
	index_inv.clear();
	total_node_num = 0;

	for (int layer_id = 0; layer_id < total_layer; ++layer_id) {
		const int source_contour_count = static_cast<int>(all_slice_points[layer_id].size());
		const bool has_containment_tree =
			layer_id < static_cast<int>(all_contour_parent_ids.size())
			&& all_contour_parent_ids[layer_id].size() == all_slice_points[layer_id].size();

		std::vector<vasco::contour_containment::ContourComponent> components;
		if (has_containment_tree) {
			contour_parent_ids[layer_id] = all_contour_parent_ids[layer_id];
			components = vasco::contour_containment::BuildMaterialComponents(
				contour_parent_ids[layer_id]);
		}
		else {
			contour_parent_ids[layer_id].assign(source_contour_count, -1);
			components.reserve(source_contour_count);
			for (int contour_id = 0; contour_id < source_contour_count; ++contour_id) {
				vasco::contour_containment::ContourComponent component;
				component.outer_contour_id = contour_id;
				components.push_back(std::move(component));
			}
		}

		const int component_count = static_cast<int>(components.size());
		total_node_num += component_count;
		slice_points[layer_id].resize(component_count);
		slice_points_contain[layer_id].resize(component_count);
		slice_points_holes[layer_id].resize(component_count);
		source_contour_ids[layer_id].resize(component_count, -1);
		source_hole_contour_ids[layer_id].resize(component_count);
		z_value[layer_id].resize(component_count);
		is_contour[layer_id].resize(component_count, false);
		adjacent_points[layer_id].resize(component_count);

		for (int component_id = 0; component_id < component_count; ++component_id) {
			const int source_id = components[component_id].outer_contour_id;
			if (source_id < 0 || source_id >= source_contour_count) {
				continue;
			}
			source_contour_ids[layer_id][component_id] = source_id;

			for (const Eigen::Vector3d& point : all_slice_points[layer_id][source_id]) {
				slice_points[layer_id][component_id].emplace_back(point.x(), point.y());
				z_value[layer_id][component_id].push_back(point.z());
			}

			if (has_containment_tree) {
				for (int hole_source_id : components[component_id].hole_contour_ids) {
					if (hole_source_id < 0 || hole_source_id >= source_contour_count) {
						continue;
					}
					source_hole_contour_ids[layer_id][component_id].push_back(hole_source_id);
					slice_points_holes[layer_id][component_id].emplace_back();
					auto& hole = slice_points_holes[layer_id][component_id].back();
					hole.reserve(all_slice_points[layer_id][hole_source_id].size());
					for (const Eigen::Vector3d& point : all_slice_points[layer_id][hole_source_id]) {
						hole.emplace_back(point.x(), point.y());
					}
				}
			}
			else if (layer_id < static_cast<int>(all_slice_points_contain.size())
				&& source_id < static_cast<int>(all_slice_points_contain[layer_id].size())
				&& !all_slice_points_contain[layer_id][source_id].empty()) {
				slice_points_holes[layer_id][component_id].emplace_back();
				auto& hole = slice_points_holes[layer_id][component_id].back();
				hole.reserve(all_slice_points_contain[layer_id][source_id].size());
				for (const Eigen::Vector3d& point : all_slice_points_contain[layer_id][source_id]) {
					hole.emplace_back(point.x(), point.y());
				}
			}

			// Preserve the old single-hole view for existing visualization code.
			if (!slice_points_holes[layer_id][component_id].empty()) {
				slice_points_contain[layer_id][component_id] =
					slice_points_holes[layer_id][component_id].front();
			}

			if (!slice_points[layer_id][component_id].empty()) {
				is_contour[layer_id][component_id] = IsContour(layer_id, component_id);
			}
			adjacent_points[layer_id][component_id].push_back(
				std::make_pair(Eigen::Vector3i(-1, -1, -1), Eigen::Vector3i(-1, -1, -1)));
			adjacent_points[layer_id][component_id].push_back(
				std::make_pair(Eigen::Vector3i(-1, -1, -1), Eigen::Vector3i(-1, -1, -1)));
		}
	}
	SetIndexMapping();
}

void Data::SetIndexMapping()
{
	int num = 0;
	for (int i = 0; i < static_cast<int>(slice_points.size()); ++i) {
		for (int j = 0; j < static_cast<int>(slice_points[i].size()); ++j) {
			index[num] = std::make_pair(i, j);
			index_inv[std::make_pair(i, j)] = num;
			++num;
		}
	}
}

bool Data::IsContour(int i, int j)
{
	if (i < 0 || i >= static_cast<int>(slice_points.size())
		|| j < 0 || j >= static_cast<int>(slice_points[i].size())
		|| slice_points[i][j].size() < 2) {
		return false;
	}

	const Eigen::Vector2d& first = slice_points[i][j].front();
	const Eigen::Vector2d& last = slice_points[i][j].back();
	return std::abs(first.x() - last.x()) <= 2
		&& std::abs(first.y() - last.y()) <= 2;
}
