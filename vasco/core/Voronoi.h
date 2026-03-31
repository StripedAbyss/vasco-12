#pragma once
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "../../vectornd.h"
#include "../../visual.h"

namespace vasco
{
    struct VoronoiCell
    {
        bool is_available{false};
        int site{-1};
        std::vector<int> adjacent_cells;
        std::vector<Eigen::Vector3d> all_points_in_polygon;
        std::vector<Eigen::Vector3d> all_normal_in_all_ori;
    };

    // 生成 Voronoi 单元
    // thresholdZ: z 高度阈值
    // visualize: 是否生成可视化 obj
    // fileName: 基础文件名（不含扩展名）
    void BuildVoronoiCells(const Eigen::MatrixXd& V,
                           const Eigen::MatrixXi& F,
                           double thresholdZ,
                           std::vector<VoronoiCell>& outCells,
                           std::vector<Eigen::Vector3d>& outBottomVertices,
                           bool visualize,
                           const std::string& fileName);
}