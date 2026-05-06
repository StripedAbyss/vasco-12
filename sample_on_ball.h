#pragma once
#ifndef SAMPLE_ON_BALL_
#define SAMPLE_ON_BALL_

#include <vector>
#include <cmath>

#include "helpers.h"
class SAMPLE_ON_BALL
{
public:
    void OrientationSamplePoints();
    void OrientationSamplePoints_2();
	std::vector<Eigen::Vector3d> sample_points; //存储球面采样点
};


#endif 