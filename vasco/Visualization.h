#pragma once
#ifndef VASCO_VISUALIZATION_H_
#define VASCO_VISUALIZATION_H_

#include <vector>
#include <string>
#include <Eigen/Core>
#include "generateClosedTube.h" // 提供 generateClosedTube 所需类型
#include "../PolygonIntersection.h"
#include "core/GeometryUtils.h"
#include "core/ObjIO.h"

namespace vasco
{
    // 生成红色球体（不可达点可视化）
    void createRedBalls(const std::string& file_name,
                        const std::vector<Eigen::MatrixXd>& vis_points);

    // 生成绿色球体（覆盖点可视化，带颜色映射）
    void createGreenBalls(const std::string& file_name,
                          const std::vector<Eigen::MatrixXd>& vis_points,
                          const std::vector<double>& color_map);

    // 从 HybridManufacturing 中迁出的分层可视化（细圆管）
    void visualize_layers(const std::string& file_name,
        const Eigen::Vector3f& vn,
        bool judge_continue_additive,
        int id_continue);

    // 从 HybridManufacturing 中迁出的阶梯层可视化（楼梯式）
    void visualize_layers_stair_case(const std::string& file_name,
        const std::string& vis_file_contain,
        const Eigen::Vector3f& vn,
        bool judge_continue_additive,
        int id_continue);
}

#endif