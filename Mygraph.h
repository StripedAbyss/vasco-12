#pragma once

#ifndef GRAPH_H_
#define GRAPH_H_


#include "edge.h"
#include "helpers.h"

#include       <set>
#include       <map>
#include    <vector>
#include    <string>
#include   <fstream>
#include <algorithm>

class MyGraph
{
public:
	MyGraph();
	~MyGraph();
	void AddEdge(int from, int to);
	void AddEdge_2(int from, int to);

	int total_node_num;
	std::vector<Edge> edges;
	vector<vector<int>>  G;
	vector<vector<int>>  G_2;
	vector<vector<int>>  G_3;
	std::vector<int>  in_degree;
	std::vector<int>  out_degree;
	std::vector<bool> node_visited, edge_deleted;
	void UpdateDegree(int node_id,int change);
	void UpdateDegree_2(int node_id, int change);
};
#endif