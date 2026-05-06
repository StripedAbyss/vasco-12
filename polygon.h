#pragma once
#ifndef POLYGON_H
#define POLYGON_H


#include <vector>

#include "helpers.h"

class Polygon
{
public:
	Polygon(const std::vector<Eigen::Vector2d>& 
	
	
	);
	~Polygon();

	bool JudgePointInside(Eigen::Vector2d p);
	void Draw(int index);
private:
	std::vector<Eigen::Vector2d> points;
	double minx, maxx, miny, maxy;
	int find = 0;
};


#endif // !POLYGON_H
