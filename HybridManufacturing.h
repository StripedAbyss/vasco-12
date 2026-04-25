#pragma once
#ifndef HYBRIDMANUFACTURING_H_
#define HYBRIDMANUFACTURING_H_
#include "polygon.h"
#include "Mygraph.h"
#include "data.h"
#include "sample_on_ball.h"
#include     <set>
#include  <vector>
#include  <string>
#include "visual.h"
#include<igl/readOBJ.h>
#include<igl/readSTL.h>
#include<igl/writeOBJ.h>
#include<igl/writeSTL.h>
#include<queue>
#include <algorithm>
#include <Eigen/Dense>
#include<iostream>
#include <cmath>
#include <regex>
#include <array>
#include <string>
#include<sstream>
#include<map>
#include <cmath>
#include "katana/katana.h"
#include"layer_graph.h"
#include<list>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>
//#include <CGAL/boost/graph/IO/OBJ.h>
//#include <CGAL/boost/graph/IO/STL.h>
#include <vector>
#include <numeric>
#include"PolygonIntersection.h"
#include <cstdio>
#include <cstdlib>

#include "vectornd.h"
#include "geometry.h"
#include "importstl.h"
#include "exportobj.h"
#include <tuple>
#include <ranges>
//#include"Point3f.h"

//#include <CGAL/intersections.h>
//#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/draw_polygon_2.h>
#include <iostream>
#include"earcut.hpp"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace std;

#include "vasco/core/Constants.h"
#include "vasco/core/Types.h"
#include "vasco/core/GeometryUtils.h"
using vasco::check_inside_2;
using vasco::core::distance3d;
using vasco::core::distance2d;
using vasco::core::distanceVec3;
using vasco::core::faceNormal;
using vasco::core::triangleArea;
using vasco::core::isAnticlockwise;
using vasco::core::checkHaveNoOtherPoint;
using vasco::core::pointsAreAnticlockwise;
using vasco::core::newGenerateMesh;

//using VEctor = std::array<double, 3>;
//using TRiangle = std::array<int, 3>;

using VEctor = vasco::core::Vec3;
using TRiangle = vasco::core::Tri3;
using area_S = vasco::core::OrientationCollision;
//using all_value = vasco::core::EvaluationScores;

#include "vasco/core/Voronoi.h"  // 新增
using vasco::VoronoiCell;        // 兼容旧代码

#include "vasco/Visualization.h"
using vasco::createRedBalls;
using vasco::createGreenBalls;
using vasco::visualize_layers;
using vasco::visualize_layers_stair_case;


#include "vasco/Slicer.h"
using  Slicer_2 = vasco::Slicer;
using  all_value = vasco::core::all_value;

#include "vasco/core/VascoUtils.h"
using vasco::SortEdges;
using vasco::quick_sort;
using vasco::Anticlockwise;
using vasco::calculate_fragile_value;
using vasco::calculate_projected_area;
using vasco::FindAllCutLayers;
using vasco::sort_candidate_nodes;

#include "vasco/core/GraphUtils.h"
using vasco::GeneralGraph_DArraySArraySpatVarying;
using vasco::GeneralGraph_DArraySArraySpatVarying2;

#include "vasco/generateClosedTube.h"

class HybridManufacturing
{
public:
	HybridManufacturing(std::string file_name, std::string suf, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd Normals);
	~HybridManufacturing();

	//void GetVoronoiCells();              // 原函数保留
    //void GetVoronoiCells1();          // 移除声明
    void InitializeVoronoi(const std::vector<VoronoiCell>& cells,
                           const std::vector<Eigen::Vector3d>& bottom); // 新增
	
	void InitializePolyscope();

	int CollisionDetectionForSubtractiveManufacturing(cutter the_nozzle);
	void GetALLFragileVertex(SAMPLE_ON_BALL sampling);
	void detect_collision_with_printing_platform(int& index,vector<int>& candidate_nodes, vector<all_value>& all_calculated_value, vector<vector<cv::Point3d>> all_cut_layers, Eigen::Vector3d ori_now, nozzle the_nozzle);
	all_value GainMesh(Slicer_2& slicer, vector<vector<cv::Point3d>> all_cut_layers, Eigen::Vector3d vector_after, int height_of_beam_search, int cont_number_of_queue, int index_of_pre_node, vector<int> all_cut_layers_dependency_layer, bool flag_is_continue_block, int id_continue);
	void CutMesh(vector<vector<cv::Point3d>> all_layers, vector<vector<cv::Point3d>> all_layers_contain,vector<vector<cv::Point3d>> all_cut_layers, Eigen::Vector3d vector_after, int height_of_beam_search, int cont_number_of_queue, int index_of_pre_node, vector<int> all_cut_layers_dependency_layer, bool& jud_outer_beam_search_terminate, vector<TRiangle>& current_remove_triangles, Slicer_2& current_slicer, bool judge_continue_additive,bool flag_is_continue_block,int pre_cont_number_of_queue, vector<bool>& jud_error, int id_node,int id_continue,vector<int> flag_cut_layers_is_hole);
	void subtractive_accessibility_decomposition(vector<TRiangle> need_detect_triangle, int height_of_beam_search, int cont_number_of_queue, cutter cutting_tool, Slicer_2 current_slicer);
	//void subtractive_accessibility_decomposition(vector<TRiangle> need_detect_triangle, int height_of_beam_search, int cont_number_of_queue, cutter cutting_tool, Slicer_2 current_slicer);
	void subtractive_accessibility_decomposition_within_2_blocks(int height_of_beam_search, cutter cutting_tool);
	vector<vector<int>> getAccessOri(const Slicer_2& slicer, Slicer_2& slicer_load_patch, vector<vasco::core::Vec3>& all_sample_points_in_triangles, cutter cutting_tool);
	void outer_beam_search(nozzle the_nozzle, cutter cutting_tool);
	void DFS_search(Layer_Graph layer_graph, bool& flag_continue, bool previous_is_continue, vector<bool> judge_S_be_searched, vector<bool> judge_covering_points_be_searched, bool& jud_admit);
	//vector<vector<vector<cv::Point3d>>> DFS_search(Layer_Graph layer_graph, vector<vector<int>>& final_pathes_include_S, vector<vector<int>>& final_pathes_include_sample_points, vector<bool> judge_S_be_searched, vector<vector<int>>& all_cut_layers, vector<vector<int>>& all_cut_layers_dependency_layer, vector<vector<int>>& final_paths_include_covering_points, vector<bool> judge_covering_points_be_searched, vector<vector<area_S>> ori_all_the_covering_points);
	void sort_candidate_nodes(vector<int>& candidate_nodes, vector<vector<vector<cv::Point3d>>> Tree_nodes, vector<vector<int>> final_pathes_include_S, vector<all_value>& all_calculated_value, vector<vector<int>> Tree_nodes_cut_layers, int pre_tree_nodes[], vector<double> Tree_nodes_larger_base, vector<vector<int>> final_pathes_include_covering_points, int height_of_beam_search, vector<vector<Eigen::Vector3d>> save_ori, vector<all_value>& pure_value, int id_continue);
	void sort_candidate_nodes(vector<int>& candidate_nodes, vector<vector<int>> Tree_nodes_for_S);
	void subtractive_remove_output(const vector<TRiangle> &need_detect_triangle, const Slicer_2 &current_slicer, int height_of_beam_search);
	bool open_vis_voronoi;
	bool open_vis_red_points;
	bool open_vis_green_points;
	bool open_change_orientation;
	bool open_vis_stair_case;

	/**
	 * @brief 检测刀具是否与指定的三角形面片发生碰撞
	 * @param center_point 刀尖球心坐标
	 * @param target_cell_vertices 目标单元的边界顶点（已旋转）
	 * @param max_z_target 目标单元边界顶点的最大z值
	 * @param nozzle 刀具参数
	 * @return true 表示发生碰撞，false 表示无碰撞
	 */
	bool CheckToolCollisionWithCell(
		const Eigen::Vector3d& center_point,
		const std::vector<Eigen::MatrixXd>& target_cell_vertices,
		double max_z_target,
		const cutter& nozzle,
		double z_threshold_divisor = 1.0,
		double xy_tolerance = 5.0) const;

	void SortCutLayersByHeight(
		vector<vector<cv::Point3d>>& all_cut_layers,
		vector<int>& all_cut_layers_dependency_layer,
		vector<int>& flag_cut_layers_is_hole,
		map<int, int>& follow_index) const;

	Slicer_2 LoadSlicerForCutMesh(
		bool flag_is_continue_block,
		int height_of_beam_search,
		int index_of_pre_node,
		int id_continue) const;

	void RotateSlicerPositions(
		Slicer_2& slicer,
		const Eigen::Vector3d& vector_before,
		const Eigen::Vector3d& vector_after) const;

	void RotateLayersForVisualization(
		vector<vector<cv::Point3d>>& all_layers,
		vector<vector<cv::Point3d>>& all_layers_contain,
		const Eigen::Vector3d& vector_after,
		const Eigen::Vector3d& vector_before) const;

	void VisualizeCutLayers(
		const vector<vector<cv::Point3d>>& all_layers,
		const vector<vector<cv::Point3d>>& all_layers_contain,
		int height_of_beam_search,
		int cont_number_of_queue,
		int index_of_pre_node,
		bool judge_continue_additive,
		int id_continue,
		const Eigen::Vector3d& vector_after) const;

private:
	Eigen::MatrixXd V;
	std::vector<Eigen::MatrixXd> V_2;
	Eigen::MatrixXi F;
	Eigen::MatrixXd Normals;
	std::vector<Eigen::Vector3d> V_bottom;              // 仍使用
	std::vector<VoronoiCell> all_voronoi_cells; // 类型改为命名空间版本

	vector<vector<bool>> is_fragile_V;	//is_fragile_V存储在各个方向下每个顶点是否为fragile点
	vector<vector<area_S>> all_the_area_S;	//all_the_area_S[i][j]存储第i个不可达点的第j个area_S，area_S包含碰撞单元和方向
	vector<vector<area_S>> all_the_covering_points;	//all_the_covering_points存储每个覆盖点对应的碰撞单元和方向area_S
	unordered_map<int, int> map_S_and_vertex;	//map_S_and_vertex存储area_S.id_to_point到原始网格顶点的映射,记录不可达点在voronoi图中的索引 <不可达点id,voronoi的id>
	unordered_map<int, int> map_S_and_vertex_inv;	//map_S_and_vertex的反向映射，记录原始voronoi单元对应的不可达点索引 <voronoi的id,不可达点id>
	unordered_map<int, int> map_covering_points_and_vertex;
	unordered_map<int, int> map_covering_points_and_vertex_inv;
	vector<int> flag_voronoi_has_been_printed;
	vector<saved_mesh> all_saved_mesh;
	SAMPLE_ON_BALL sampling_subtractive;
	//SamplePoints AllSamplePoints;
	//vector<VoronoiCell> all_voronoi_cells;

	vector<vector<int>> ori_num_points_of_ori_in_all_the_area_S;
	vector<vector<int>> pathes_include_S, pathes_include_sample_points, paths_include_covering_points;
	vector<vector<int>> all_cut_layers; // 存储每条路径中的切割层在final_pathes[i]中的下标
	vector<vector<int>> all_cut_layers_dependency_layer; // 存储每条路径中每个切割层所依赖的切割层数量
	vector<vector<area_S>> ori_all_the_covering_points;
	vector<vector<vector<cv::Point3d>>> all_solutions_of_selected_layers;
	vector<vector<vector<cv::Point3d>>> all_solutions_of_selected_layers_contain;

	double time_build_subtractive_graph;
	int num_inaccessible_points;
	int cont_extra_additive_orientation;
	string file_name;
	string suf;

	std::vector<std::vector<int>> EvaluateMergedPatchToolCollision(
		const Slicer_2& merged_patch,
		const std::vector<int>& merged_face_source_patch_id,
		cutter cutting_tool) const;

	Slicer_2 MergeBlockPatchesWithDedup(
		int max_patch_index,
		std::vector<int>& merged_vertex_source_patch_id,
		std::vector<int>& merged_face_source_patch_id,
		double merge_eps = 1e-6) const;

	void ExportMergedPatchFaceColorOBJ(
		const Slicer_2& merged_patch,
		const std::vector<int>& max_collision_patch_per_face,
		const std::string& color_obj_file) const;
};




#endif