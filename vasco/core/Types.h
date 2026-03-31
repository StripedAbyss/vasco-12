#pragma once
#include <array>
#include <vector>

namespace vasco::core {

// 基础向量与三角面索引
using Vec3    = std::array<double, 3>;
using Tri3    = std::array<int, 3>;

// 原 area_S => 更语义化命名
struct OrientationCollision {
    int pointId;   // 原 id_to_point
    int oriId;     // 原 id_ori
    OrientationCollision(int p = 0, int o = 0) : pointId(p), oriId(o) {}
};

// 原 all_value => 评估指标集合
struct EvaluationScores {
    double baseLarge              = 0.0;
    double selfSupport            = 0.0;
    double areaS                  = 0.0;
    double sliceLayers            = 0.0;
    double clippingPlanes         = 0.0;
    double distanceToFragile      = 0.0;
    double coveringPoints         = 0.0;
    double fragile                = 0.0;
    double remainingFaces         = 0.0;
    double orientationSimilarity  = 0.0;
    double projectedRatio         = 0.0;
};

struct all_value
{
    double large_base;
    double value_of_self_support;
    double value_of_area_S;
    double value_of_more_slice_layers;
    double value_of_less_clipping_plane;
    double value_of_distance_to_fragile_regions;
    double value_of_covering_points;
    double value_of_fragile;
    double number_of_remaining_face;
    double value_of_orientation;
    double value_of_projected;
};

const double MAX_D = 1e18;

const int NUM_COLORS_HZ = 21;

const Tri3 COLORS_HZ[21] = {
    {223,  70,  97}, {241, 128,  58}, {255, 205,   0}, {108, 194,  74}, {  0, 174, 199}, {  0, 114, 206}, {135,  24, 157},
    {172,  20,  90}, {190,  77,   0}, {218, 170,   0}, {  0, 122,  62}, {  0, 140, 149}, {  0,  71, 187}, {117,  59, 189},
    {246, 117, 153}, {255, 178,  91}, {239, 223,   0}, {151, 215,   0}, {  5, 195, 221}, {108, 172, 228}, {172,  79, 198}
};

const Vec3 COLORS_HZ_FLOAT[21] = {
    {223.0 / 255.0,  70.0 / 255.0,  97.0 / 255.0}, {241.0 / 255.0, 128.0 / 255.0,  58.0 / 255.0}, {255.0 / 255.0, 205.0 / 255.0,   0.0 / 255.0},
    {108.0 / 255.0, 194.0 / 255.0,  74.0 / 255.0}, {  0.0 / 255.0, 174.0 / 255.0, 199.0 / 255.0}, {  0.0 / 255.0, 114.0 / 255.0, 206.0 / 255.0},
    {135.0 / 255.0,  24.0 / 255.0, 157.0 / 255.0}, {172.0 / 255.0,  20.0 / 255.0,  90.0 / 255.0}, {190.0 / 255.0,  77.0 / 255.0,   0.0 / 255.0},
    {218.0 / 255.0, 170.0 / 255.0,   0.0 / 255.0}, {  0.0 / 255.0, 122.0 / 255.0,  62.0 / 255.0}, {  0.0 / 255.0, 140.0 / 255.0, 149.0 / 255.0},
    {  0.0 / 255.0,  71.0 / 255.0, 187.0 / 255.0}, {117.0 / 255.0,  59.0 / 255.0, 189.0 / 255.0}, {246.0 / 255.0, 117.0 / 255.0, 153.0 / 255.0},
    {255.0 / 255.0, 178.0 / 255.0,  91.0 / 255.0}, {239.0 / 255.0, 223.0 / 255.0,   0.0 / 255.0}, {151.0 / 255.0, 215.0 / 255.0,   0.0 / 255.0},
    {  5.0 / 255.0, 195.0 / 255.0, 221.0 / 255.0}, {108.0 / 255.0, 172.0 / 255.0, 228.0 / 255.0}, {172.0 / 255.0,  79.0 / 255.0, 198.0 / 255.0}
};

} // namespace vasco::core