#include "data.h"

void Data::ReadData(vector<vector<vector<Vertex>>> all_slice_points, vector<vector<vector<Vertex>>> all_slice_points_contain)
{
	int total_layer, s_num, p_num;
	double x, y, z, x2, y2;
	total_layer = all_slice_points.size();
	this->slice_points.resize(total_layer);
	this->slice_points_contain.resize(total_layer);
	this->z_value.resize(total_layer);
	this->is_contour.resize(total_layer);
	this->adjacent_points.resize(total_layer);
	this->total_node_num = 0;
	for (int i = 0; i < total_layer; i++) {
		s_num = all_slice_points[i].size();
		this->total_node_num += s_num;
		this->slice_points[i].resize(s_num);
		this->slice_points_contain[i].resize(s_num);
		this->z_value[i].resize(s_num);
		this->is_contour[i].resize(s_num);
		this->adjacent_points[i].resize(s_num);
		for (int j = 0; j < s_num; j++) {
			p_num = all_slice_points[i][j].size();
			for (int k = 0; k < p_num; k++) {
				x = all_slice_points[i][j][k].x;
				y = all_slice_points[i][j][k].y;
				z = all_slice_points[i][j][k].z;
				this->slice_points[i][j].push_back(Eigen::Vector2d(x, y));
				this->z_value[i][j].push_back(z);
			}
			if (all_slice_points_contain[i][j].size() != 0) {
				for (int k = 0; k < all_slice_points_contain[i][j].size(); k++) {
					x2 = all_slice_points_contain[i][j][k].x;
					y2 = all_slice_points_contain[i][j][k].y;
					this->slice_points_contain[i][j].push_back(Eigen::Vector2d(x2, y2));
				}
			}

			// avoid the segment only has one point
			if (p_num == 1) {
				this->slice_points[i][j].push_back(Eigen::Vector2d(x + 1, y + 1));
				this->z_value[i][j].push_back(z);
			}
			if (IsContour(i, j)) {

				this->is_contour[i][j] = true;
			}
			this->adjacent_points[i][j].push_back(std::make_pair(Eigen::Vector3i(-1, -1, -1), Eigen::Vector3i(-1, -1, -1)));
			this->adjacent_points[i][j].push_back(std::make_pair(Eigen::Vector3i(-1, -1, -1), Eigen::Vector3i(-1, -1, -1)));
		}
	}
	SetIndexMapping();
}

void Data::SetIndexMapping()
{
	int num = 0;
	for (int i = 0; i < slice_points.size(); i++) {
		for (int j = 0; j < slice_points[i].size(); j++) {
			index[num] = std::make_pair(i, j);
			index_inv[std::make_pair(i, j)] = num;
			num++;
		}
	}
}

bool Data::IsContour(int i, int j)
{

	Eigen::Vector2d p1 = this->slice_points[i][j][0];
	Eigen::Vector2d p2 = this->slice_points[i][j][this->slice_points[i][j].size() - 1];
	//std::cout << i << " " << j << " " << p1 << " " << p2 << std::endl;
	//if (JudgePointEqual(p1, p2)) return true;
	if (std::abs(p1.x() - p2.x()) <= 2 && std::abs(p1.y() - p2.y()) <= 2) return true;
	else return false;
}
