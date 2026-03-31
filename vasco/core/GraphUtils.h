#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
// 从 vasco/core 到工程根目录的 gco-v3.0
#include "../../gco-v3.0/GCoptimization.h"


namespace vasco {

    // 通用图割（General Graph Cut）
    // - num_pixels: 像素/节点数量
    // - num_labels: 标签数量
    // - data_value: 数据项代价，data_value[i][l] 为像素 i 的标签 l 的代价
    // - pixels_relations: 邻接边列表，形如 {{u, v}, ...}
    // - length_edges: 与 pixels_relations 一一对应的边权
    std::vector<int> GeneralGraph_DArraySArraySpatVarying(
        int num_pixels,
        int num_labels,
        const std::vector<std::vector<int>>& data_value,
        const std::vector<std::vector<int>>& pixels_relations,
        const std::vector<int>& length_edges) noexcept;

    std::vector<int> GeneralGraph_DArraySArraySpatVarying2(
        int num_pixels, 
        int num_labels, 
        const std::vector<std::vector<int>>& data_value, 
        const std::vector<std::pair<int, int>>& pixels_relations,
        const std::vector<int>& length_edges) noexcept;

    // 通用图割（General Graph Cut）
    // - num_pixels: 像素/节点数量
    // - num_labels: 标签数量
    // - data_value: 数据项代价，data_value[i][l] 为像素 i 的标签 l 的代价

} // namespace vasco