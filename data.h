#pragma once
#ifndef DATA_H_
#define DATA_H_


#include       <map>
#include   <fstream>

#include "helpers.h"
#include "katana/datastructures.h"
class Data
{
public:
	void ReadData(vector<vector<vector<Vertex>>> all_slice_points, vector<vector<vector<Vertex>>> all_slice_points_contain);

	int total_node_num;
	std::vector<std::vector<std::vector<Eigen::Vector2d>>> slice_points;
	std::vector<std::vector<std::vector<double>>> z_value;
	std::vector<std::vector<bool>> is_contour;
	std::map<int, std::pair<int, int>> index; //1d->2d
	std::map<std::pair<int, int>, int> index_inv; //2d->1d

	std::vector<std::vector<std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3i>>>> adjacent_points;

	std::vector<std::vector<std::vector<Eigen::Vector2d>>> slice_points_contain;

private:
	void SetIndexMapping();
	bool IsContour(int i, int j);
};

#endif