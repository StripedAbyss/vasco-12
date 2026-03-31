#pragma once

#include "Types.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
#include "../Slicer.h"
#include "GeometryUtils.h"


namespace vasco {
	// 需要拿到 area_S 的完整定义（含 oriId 字段）。
	using area_S = core::OrientationCollision;

	// 将每个子向量按元素的 oriId 升序排序。
	inline void SortEdges(std::vector<std::vector<area_S>>& edges) noexcept;

	// 计算制造需求分值：三角形平均高度越低，值越高。
	// maxZ == minZ 时返回 0 以避免除零。
	double calculateManufacturingRequire(const Slicer& slicer,
		int triangleIndex,
		double minZ,
		double maxZ) noexcept;

	// 按 min_z_triangle 的值对并行数组进行原地快速排序（与旧实现等价）。
	void quick_sort(
		std::vector<vasco::core::Tri3>& candidate_triangles,
		int left,
		int right,
		std::vector<int>& id_triangles,
		std::vector<double>& min_z_triangle,
		std::vector<vasco::core::Vec3>& min_z_point) noexcept;

	// 将 real_cutting_plane_triangles 的每个环的顶点索引按 all_slicer.positions
	// 的投影方向修正为逆时针。如判断为顺时针则反转。
	// 要求 all_slicer.positions 按 [x,y,z] 访问，且索引有效。
	void Anticlockwise(std::vector<std::vector<int>>& real_cutting_plane_triangles,
		const Slicer& all_slicer) noexcept;

	// 多边形耳切三角化：输入顶点索引序列（按 isAnti 指定的顺/逆时针），返回三角形顶点索引序列（每3个为一面）。
	std::vector<int> poufen(const Slicer& slicer,
		const std::vector<int>& index_of_points,
		bool isAnti) noexcept;

	// 统计每个切割层内的脆弱点最低距离是否在阈值内，累加计数写入 all_calculated_value.value_of_fragile，用来计算得分
	void calculate_fragile_value(
		vasco::core::all_value& all_calculated_value,
		const std::vector<std::vector<cv::Point3d>>& all_cut_layers,
		const std::vector<Eigen::MatrixXd>& Tree_nodes_fragile_V) noexcept;

	// 计算所有切割层的平均投影面积比值（候选三角形投影的凸包面积 / 切割界面面积），用来计算得分
	double calculate_projected_area(
		const Slicer& all_slicer,
		const std::vector<std::vector<vasco::core::Tri3>>& all_furcation_of_blocks,
		const std::vector<std::vector<cv::Point3d>>& all_cut_layers) noexcept;

	// 从分层依赖图路径中找出作为切割层的 layer 下标，并输出依赖层数量
	std::vector<std::vector<int>> FindAllCutLayers(
		const Layer_Graph& layer_graph,
		const std::vector<std::vector<int>>& final_pathes,
		std::vector<std::vector<int>>& all_cut_layers_dependency_layer,
		bool& jud_admit) noexcept;

	// 重载一：完整排序（评分聚合、归一化与加权）
	void sort_candidate_nodes(
		std::vector<int>& candidate_nodes,
		const std::vector<std::vector<std::vector<cv::Point3d>>>& Tree_nodes,
		const std::vector<std::vector<int>>& final_pathes_include_S,
		std::vector<vasco::core::all_value>& all_calculated_value,
		const std::vector<std::vector<int>>& Tree_nodes_cut_layers,
		const int pre_tree_nodes[],
		const std::vector<double>& Tree_nodes_larger_base,
		const std::vector<std::vector<int>>& final_pathes_include_covering_points,
		int height_of_beam_search,
		const std::vector<std::vector<Eigen::Vector3d>>& save_ori,
		std::vector<vasco::core::all_value>& pure_value,
		int id_continue) noexcept;

	// 重载二：简单排序（按 S 数量）
	void sort_candidate_nodes(
		std::vector<int>& candidate_nodes,
		const std::vector<std::vector<int>>& Tree_nodes_for_S) noexcept;



} // namespace vasco