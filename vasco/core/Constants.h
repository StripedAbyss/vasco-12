#pragma once
#include <cstddef>

namespace vasco::core {

// 数值精度统一
inline constexpr double EPS        = 1e-9;
inline constexpr double SMALL_EPS  = 1e-12;

// 常用数学常量
inline constexpr double PI         = 3.14159265358979323846;

// 全局阈值集中（后续可改为配置加载）
struct GlobalThresholds {
    double fragileDotThreshold = 0.95;     // 脆弱点法向与Z轴点积阈值
    double slicingMinPlaneGap  = 1e-4;     // 切平面高度比较容差
    double minTriangleArea     = 1e-8;     // 极小三角形过滤
};
inline GlobalThresholds THRESHOLDS;

// 搜索 / 算法参数（后续放配置文件）
struct AlgorithmParams {
    int    terminateRemainingFaces = 800;  // 示例：终止条件最少剩余面数
    double selfSupportSlopeRad     = PI/3.6;
    double ratioNotSelfSupportMax  = 0.05;
};
inline AlgorithmParams PARAMS;

} // namespace vasco::core