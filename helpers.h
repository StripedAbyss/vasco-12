#pragma once
#ifndef HELPERS_H_
#define HELPERS_H_

#include <string>
#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include <opencv2/core/core.hpp>
using namespace std;


struct nozzle
{
	double upper_surface_r;
	double lowwer_surface_r;
	double nozzle__H_total;
	double nozzle_H_half;
};

struct cutter
{
	double cylinder_r;
	double cylinder_height;
	double ball_r;
	double carriage_r;
	double carriage_height;
	double cylinder_height_threshold;
	double carriage_check_radius_sq;
	double cylinder_check_radius_sq;
	double cylinder_r_sq;
	double carriage_r_sq;
	double total_height;
};

const int terminate_threshold_of_number_of_faces = 1000;  //1500
const int    maxn = 10000;
const int    MAX_I = 100000000;
const double eps = 1e-10;
const double dependence_offset = 2;   //0.5(FDM)     //3.5(ceramic)  //5
const double MAX_D = 1e18;
const double MIN_D = -1e18;
const double dh = 2.0;
const int num_ori_sample = 200; //200

const std::string file_name = ".\\data"; //layer_graph.cpp里面用到了 引用了Mygraph

extern std::vector<std::vector<bool>> is_flatten_area;
extern std::vector<std::pair<int,int>> index_flatten_layer;
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

#endif