#pragma once
#ifndef LAYER_GRAPG_H_
#define LAYER_GRAPG_H_


#include "polygon.h"
#include "Mygraph.h"
#include "data.h"
#include "sample_on_ball.h"
#include     <set>
#include  <vector>
#include  <string>
#include "visual.h"
#include<igl/readOBJ.h>
#include<queue>
#include <algorithm>
using namespace std;

struct temp_area_S
{
	int* id_layers;
	int* id_points;
	int* id_ori;
	vector<pair<int, int>> subtractive_collision_edges;
	temp_area_S(int num) {
		id_layers = new int[num];
		id_points = new int[num];
		id_ori = new int[num];
		subtractive_collision_edges.resize(num);
	}
};

struct Area_S
{
	int id_layer;
	vector<int> id_dep_layers;
	vector<vector<pair<int, int>>> K_and_ori;
	vector<int> all_k;
};

struct point_subtractive_map
{
	int id_k;
	bool flag_ori_unaccessible[num_ori_sample];
	point_subtractive_map(int a) {
		id_k = a;
		for (int i = 0; i < num_ori_sample; i++)
			flag_ori_unaccessible[i] = false;
	}
};

class Layer_Graph : public MyGraph
{
public:
	Layer_Graph(const Data& data);
	~Layer_Graph();
	void GetTrianglesForLayers(vector<vector<vector<Vertex>>> all_slice_points, std::vector<map<pair<Vertex, Vertex>, Triangle*>> map_segment_triangles, vector<Vertex> all_vertex, Eigen::Vector3d vectorAfter, int height_of_beam_search, int id_continue);
	void GenerateDependencyEdges();
	void BuildLayerGraph(nozzle the_nozzle);
	void BuildDependencyGraph(std::vector<cv::Point3d>& all_points);
	void GetInitialOPP();
	//void MappingBackLayers(vector<Eigen::Matrix3d> all_rotMatrix); //可视化用的？
	void CollisionDetectionForAdditiveManufacturing(nozzle the_nozzle);
	void CollisionDetectionForSubtractiveManufacturing(nozzle the_nozzle);
	void creat_ball(string file_name, std::vector<Eigen::MatrixXd> vis_points);
	void MergeLayersWithDefaultOrder();
	void MergeLayersWithGreedyAlgorithm(int* cont_layers_of_patches);
	void MergeLayersWithBeamSearch(int* cont_layers_of_patches);
	void MergeLayersWithHeuristicSearch(int* cont_layers_of_patches);

	bool IsDepend_collision(int i, int j);
	void DFS(int u, std::vector<int>& initial_opp);
	void OutputInitialOpp(const std::string& file_name);
	bool compare_two_node(std::vector<int>a, std::vector<int>b);
	void DFS_One(int u, std::vector<int>& medium_path, int num_blocks);
	void BFS(std::vector<int>& medium_path, int num_blocks);
	void DFS_ALL(int fa, int u, std::vector<int>& medium_path, int num_blocks);

	Data data;
	std::vector<std::vector<int>> initial_opp_info;
	vector<pair<int, int>> temp_edges;
	//vector<pair<int, int>> temp_subtractive_edges;
	int cont_normal_dependency_edges;
	vector<vector<Triangle*>> all_triangles_of_layers;
	vector<bool> is_the_layer_self_suppot;

private:
	std::vector<std::vector<std::vector<Eigen::MatrixXd>>> current_layers;
	int* cont_layers_of_patches;
	int* cont_nodes_of_patches;
	int num_patches;
	vector<Area_S> all_the_area_S;
	int sum_layers;
	bool* flag_layer_is_accessible;
	map<int, int> map_index_layers;
	vector<Area_S> current_area_S;
	vector<vector<point_subtractive_map>> current_point_map;
	bool** has_subtractive_collision_dependency;
	vector<vector<int>> all_blocks;
	std::vector<std::vector<std::vector<int>>> all_solutions;
	
	int node;
	int* pre_tree_index;
	std::vector<std::vector<int>> tree_node;
	std::vector<std::vector<int>> temp_tree_node;
	std::queue<int> Q;
	std::queue<int> Q_2;
};


inline static cv::Point2d GetNormal(cv::Point2d p1, cv::Point2d p2) //ConstructPolygonPoints和ConstructPolygonPoints_2里面用到了
{
	cv::Point2d res;
	res.x = -(p2.y - p1.y);
	res.y = (p2.x - p1.x);
	double norm = std::sqrt(res.dot(res));
	if (norm <= eps) {
		//std::cout << "divide 0 occur !!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		//std::cout << p1 << " " << p2 << std::endl;
		return res;
	}
	res.x /= norm;
	res.y /= norm;
	return res;
}

inline static double Distance2D(cv::Point p1, cv::Point p2) { //slicer里面用到了
	return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

inline static std::vector<cv::Point2d> ConstructPolygonPoints(const std::vector<cv::Point2d>& points, double offset) { //layer_graph里面用了
	cv::Point2d dir;
	std::vector<cv::Point2d> polygon_Points;
	std::vector<cv::Point2d> inside;
	std::vector<cv::Point2d> outside;
	for (int k = 0; k < points.size(); k++) {
		if (k == points.size() - 1) {
			dir = GetNormal(points[points.size() - 2], points[points.size() - 1]);
			inside.push_back(points[points.size() - 1] + dir * offset);
		}
		else {
			dir = GetNormal(points[k], points[k + 1]);
			inside.push_back(points[k] + dir * offset);
		}
	}
	for (int k = points.size() - 1; k >= 0; k--) {
		if (k == 0) {
			dir = GetNormal(points[1], points[0]);
			outside.push_back(points[0] + dir * offset);
		}
		else {
			dir = GetNormal(points[k], points[k - 1]);
			outside.push_back(points[k] + dir * offset);
		}
	}
	polygon_Points.push_back(points[0]);
	for (int k = 0; k < inside.size(); k++) polygon_Points.push_back(inside[k]);
	polygon_Points.push_back(points[points.size() - 1]);
	for (int k = 0; k < outside.size(); k++) polygon_Points.push_back(outside[k]);
	return polygon_Points;
}


inline static std::vector<cv::Point2d> ConstructPolygonPoints_2(std::vector<cv::Point2d>& points, double offset) { //layer_graph里面用了
	cv::Point2d dir;
	std::vector<cv::Point2d> polygon_Points;
	std::vector<cv::Point2d> inside;
	std::vector<cv::Point2d> outside;
	for (int k = 0; k < points.size(); k++) {
		if (k == points.size() - 1) {
			dir = GetNormal(points[points.size() - 2], points[points.size() - 1]);
			inside.push_back(points[points.size() - 1] + dir * offset);
		}
		else {
			dir = GetNormal(points[k], points[k + 1]);
			inside.push_back(points[k] + dir * offset);
		}
	}
	for (int k = points.size() - 1; k >= 0; k--) {
		if (k == 0) {
			dir = GetNormal(points[1], points[0]);
			outside.push_back(points[0] + dir * offset);
		}
		else {
			dir = GetNormal(points[k], points[k - 1]);
			outside.push_back(points[k] + dir * offset);
		}
	}
	polygon_Points.push_back(points[0]);
	for (int k = 0; k < inside.size(); k++) polygon_Points.push_back(inside[k]);
	polygon_Points.push_back(points[points.size() - 1]);
	for (int k = 0; k < outside.size(); k++) polygon_Points.push_back(outside[k]);
	return polygon_Points;
}

#endif // !LAYER_GRAPG_H_