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
	void MappingBackLayers(vector<Eigen::Matrix3d> all_rotMatrix);
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

#endif // !LAYER_GRAPG_H_