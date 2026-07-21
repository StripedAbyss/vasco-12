#pragma once
#ifndef DATA_H_
#define DATA_H_


#include       <map>
#include   <fstream>

#include "helpers.h"
#include "katana/datastructures.h"

using std::vector;

class Data
{
public:
	void ReadData(
		vector<vector<vector<Eigen::Vector3d>>> all_slice_points,
		vector<vector<vector<Eigen::Vector3d>>> all_slice_points_contain,
		const vector<vector<int>>& all_contour_parent_ids = {});

	int total_node_num;
	std::vector<std::vector<std::vector<Eigen::Vector2d>>> slice_points;
	std::vector<std::vector<std::vector<double>>> z_value;
	std::vector<std::vector<bool>> is_contour;
	std::map<int, std::pair<int, int>> index; //1d->2d
	std::map<std::pair<int, int>, int> index_inv; //2d->1d

	std::vector<std::vector<std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3i>>>> adjacent_points;

	// Compatibility view used by old single-hole visualization code.
	std::vector<std::vector<std::vector<Eigen::Vector2d>>> slice_points_contain;
	// All direct holes belonging to each even-depth material component.
	std::vector<std::vector<std::vector<std::vector<Eigen::Vector2d>>>> slice_points_holes;
	// Maps graph components back to raw slice-contour IDs.
	std::vector<std::vector<int>> source_contour_ids;
	std::vector<std::vector<std::vector<int>>> source_hole_contour_ids;
	// Parent-only containment trees in the raw contour index space.
	std::vector<std::vector<int>> contour_parent_ids;

private:
	void SetIndexMapping();
	bool IsContour(int i, int j);
};

#endif
