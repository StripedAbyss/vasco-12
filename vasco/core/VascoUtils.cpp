#include "VascoUtils.h"
#include <cassert>

namespace vasco {

    void SortEdges(std::vector<std::vector<area_S>>& edges) noexcept
    {
        for (auto& group : edges)
        {
            std::sort(group.begin(), group.end(),
                [](const area_S& a, const area_S& b) { return a.oriId < b.oriId; });
        }
    }

    double calculateManufacturingRequire(const Slicer& slicer,
        int triangleIndex,
        double minZ,
        double maxZ) noexcept
    {
        if (triangleIndex < 0 || triangleIndex >= static_cast<int>(slicer.triangles.size())) {
            std::cout << "if (triangleIndex < 0 || triangleIndex >= static_cast<int>(slicer.triangles.size())) {";
            return 0.0;       
        }
        if (maxZ - minZ == 0.0) {
			std::cout << "calculateManufacturingRequire (maxZ - minZ == 0.0) {";
            assert(0);
            return 0.0;
        }

        const auto& tri = slicer.triangles[triangleIndex];
        double z0 = slicer.positions[tri[0]][2];
        double z1 = slicer.positions[tri[1]][2];
        double z2 = slicer.positions[tri[2]][2];
        double averageZ = (z0 + z1 + z2) / 3.0;
        return (maxZ - averageZ) / (maxZ - minZ);
    }


    void quick_sort(std::vector<core::Tri3>& candidate_triangles,
        int left,
        int right,
        std::vector<int>& id_triangles,
        std::vector<double>& min_z_triangle,
        std::vector<core::Vec3>& min_z_point) noexcept
    {
        if (left > right) {
            return;
        }

        const double pivotZ = min_z_triangle[left];
        auto pivotTri = candidate_triangles[left];
        const int pivotId = id_triangles[left];
        auto pivotMinPt = min_z_point[left];

        int i = left;
        int j = right;

        while (i != j) {
            while (i < j && min_z_triangle[j] >= pivotZ) {
                --j;
            }
            while (i < j && min_z_triangle[i] <= pivotZ) {
                ++i;
            }
            if (i < j) {
                std::swap(candidate_triangles[i], candidate_triangles[j]);
                std::swap(id_triangles[i], id_triangles[j]);
                std::swap(min_z_triangle[i], min_z_triangle[j]);
                std::swap(min_z_point[i], min_z_point[j]);
            }
        }

        min_z_triangle[left] = min_z_triangle[i];
        min_z_triangle[i] = pivotZ;
        candidate_triangles[left] = candidate_triangles[i];
        candidate_triangles[i] = pivotTri;
        id_triangles[left] = id_triangles[i];
        id_triangles[i] = pivotId;
        min_z_point[left] = min_z_point[i];
        min_z_point[i] = pivotMinPt;

        quick_sort(candidate_triangles, left, i - 1, id_triangles, min_z_triangle, min_z_point);
        quick_sort(candidate_triangles, i + 1, right, id_triangles, min_z_triangle, min_z_point);
    }

    void Anticlockwise(std::vector<std::vector<int>>& real_cutting_plane_triangles, 
        const Slicer& all_slicer) noexcept
    {
        for (auto& ring : real_cutting_plane_triangles) {
            const int n = static_cast<int>(ring.size());
            if (n < 3) {
                continue;
            }

            // 选取 x 最大的点作为参考点
            double max_x = -std::numeric_limits<double>::infinity();
            int idx = 0;
            for (int i = 0; i < n; ++i) {
                const auto vi = ring[i];
                const double x = all_slicer.positions[vi][0];
                if (x > max_x) {
                    max_x = x;
                    idx = i;
                }
            }

            auto at = [&](int i) -> const std::array<double, 3>&{
                return all_slicer.positions[ring[(i + n) % n]];
                };

            // 以 x 最大点为中心做二维叉积，d>0 认为顺时针，需反转为逆时针
            const auto& p = all_slicer.positions[ring[idx]];
            const auto& pL = at(idx - 1);
            const auto& pR = at(idx + 1);

            const double d = (p[0] - pL[0]) * (p[1] - pR[1])
                           - (p[1] - pL[1]) * (p[0] - pR[0]);
            if (d > 0.0) {
                std::reverse(ring.begin(), ring.end());
            }
        }
    }

    std::vector<int> poufen(const Slicer& slicer,
        const std::vector<int>& index_of_points,
        bool isAnti) noexcept
    {

        std::vector<int> returnIndex;
        if (index_of_points.size() < 3) {
            return returnIndex;
        }

        // 构建点列表（按 isAnti 控制方向）
        std::list<cv::Point3f> pointList;
        std::list<int> pointIndex;

        if (!isAnti) {
            for (int i = static_cast<int>(index_of_points.size()) - 1; i >= 0; --i) {
                const int idx = index_of_points[static_cast<size_t>(i)];
                cv::Point3f temp_point(slicer.positions[idx][0], slicer.positions[idx][1], slicer.positions[idx][2]);
				pointList.push_back(temp_point);
                //pointList.push_back({ slicer.positions[idx][0], slicer.positions[idx][1], slicer.positions[idx][2] });
                pointIndex.push_back(idx);
            }
        }
        else {
            for (int i = 0; i < static_cast<int>(index_of_points.size()); ++i) {
                const int idx = index_of_points[static_cast<size_t>(i)];
				cv::Point3f temp_point(slicer.positions[idx][0], slicer.positions[idx][1], slicer.positions[idx][2]);
				pointList.push_back(temp_point);
                //pointList.push_back({ slicer.positions[idx][0], slicer.positions[idx][1], slicer.positions[idx][2] });
                pointIndex.push_back(idx);
            }
        }

        if (pointList.size() < 3) {
            return returnIndex;
        }

        auto iter = pointList.begin();
        auto iter_2 = pointIndex.begin();

        auto next_wrap = [&](auto it, int n = 1) {
            auto tmp = it;
            for (int k = 0; k < n; ++k) {
                tmp = std::next(tmp);
                if (tmp == pointList.end()) tmp = pointList.begin();
            }
            return tmp;
            };
        auto next_wrap_i = [&](auto it, int n = 1) {
            auto tmp = it;
            for (int k = 0; k < n; ++k) {
                tmp = std::next(tmp);
                if (tmp == pointIndex.end()) tmp = pointIndex.begin();
            }
            return tmp;
            };

        auto p1 = *iter;
        auto p2 = *next_wrap(iter);
        auto p3 = *next_wrap(iter, 2);
        int t1 = *iter_2;
        int t2 = *next_wrap_i(iter_2);
        int t3 = *next_wrap_i(iter_2, 2);

        while (pointList.size() >= 3) {
            int ccw = core::isAnticlockwise(p1, p2, p3); // 1: 逆时针；2: 共线；0: 顺时针
            if (ccw == 1 && core::checkHaveNoOtherPoint(p1, p2, p3, pointList)) {
                // 输出一个三角面
                returnIndex.push_back(t1);
                returnIndex.push_back(t2);
                returnIndex.push_back(t3);

                // 删除耳朵的中点 p2
                auto delP = std::next(iter);
                auto delI = std::next(iter_2);
                if (delP == pointList.end()) {
                    delP = pointList.begin();
                    delI = pointIndex.begin();
                }
                pointList.erase(delP);
                pointIndex.erase(delI);

                if (pointList.size() >= 3) {
                    if (std::next(iter) == pointList.end()) {
                        p2 = *pointList.begin();
                        p3 = *std::next(pointList.begin());
                        t2 = *pointIndex.begin();
                        t3 = *std::next(pointIndex.begin());
                    }
                    else if (std::next(iter, 2) == pointList.end()) {
                        p2 = *std::next(iter);
                        p3 = *pointList.begin();
                        t2 = *std::next(iter_2);
                        t3 = *pointIndex.begin();
                    }
                    else {
                        p2 = *std::next(iter);
                        p3 = *std::next(iter, 2);
                        t2 = *std::next(iter_2);
                        t3 = *std::next(iter_2, 2);
                    }
                }
            }
            else if (ccw == 2) {
                // 共线，删除中点 p2
                auto delP = std::next(iter);
                auto delI = std::next(iter_2);
                if (delP == pointList.end()) {
                    delP = pointList.begin();
                    delI = pointIndex.begin();
                }
                pointList.erase(delP);
                pointIndex.erase(delI);

                if (pointList.size() >= 3) {
                    if (std::next(iter) == pointList.end()) {
                        p2 = *pointList.begin();
                        p3 = *std::next(pointList.begin());
                        t2 = *pointIndex.begin();
                        t3 = *std::next(pointIndex.begin());
                    }
                    else if (std::next(iter, 2) == pointList.end()) {
                        p2 = *std::next(iter);
                        p3 = *pointList.begin();
                        t2 = *std::next(iter_2);
                        t3 = *pointIndex.begin();
                    }
                    else {
                        p2 = *std::next(iter);
                        p3 = *std::next(iter, 2);
                        t2 = *std::next(iter_2);
                        t3 = *std::next(iter_2, 2);
                    }
                }
            }
            else {
                // 顺时针，窗口前移
                ++iter;
                ++iter_2;
                if (iter == pointList.end()) {
                    iter = pointList.begin();
                    iter_2 = pointIndex.begin();
                    p1 = *iter;
                    p2 = *next_wrap(iter);
                    p3 = *next_wrap(iter, 2);
                    t1 = *iter_2;
                    t2 = *next_wrap_i(iter_2);
                    t3 = *next_wrap_i(iter_2, 2);
                }
                else if (std::next(iter) == pointList.end()) {
                    p1 = *iter;
                    p2 = *pointList.begin();
                    p3 = *std::next(pointList.begin());
                    t1 = *iter_2;
                    t2 = *pointIndex.begin();
                    t3 = *std::next(pointIndex.begin());
                }
                else if (std::next(iter, 2) == pointList.end()) {
                    p1 = *iter;
                    p2 = *std::next(iter);
                    p3 = *pointList.begin();
                    t1 = *iter_2;
                    t2 = *std::next(iter_2);
                    t3 = *pointIndex.begin();
                }
                else {
                    p1 = *iter;
                    p2 = *std::next(iter);
                    p3 = *std::next(iter, 2);
                    t1 = *iter_2;
                    t2 = *std::next(iter_2);
                    t3 = *std::next(iter_2, 2);
                }
            }
        }

        return returnIndex;
    }

    void calculate_fragile_value(
        vasco::core::all_value& all_calculated_value,
        const std::vector<std::vector<cv::Point3d>>& all_cut_layers,
        const std::vector<Eigen::MatrixXd>& Tree_nodes_fragile_V) noexcept
    {
        constexpr double threshold = 6.0;
        double sum_value_of_fragile = 0.0;

        for (size_t t = 0; t < all_cut_layers.size(); ++t) {
            if (all_cut_layers[t].empty()) {
                continue;
            }
            double boundary_bottom = 999999.0;
            double boundary_left = 999999.0;
            double boundary_top = -999999.0;
            double boundary_right = -999999.0;

            for (const auto& p : all_cut_layers[t]) {
                boundary_top = std::max(boundary_top, p.y);
                boundary_bottom = std::min(boundary_bottom, p.y);
                boundary_right = std::max(boundary_right, p.x);
                boundary_left = std::min(boundary_left, p.x);
            }

            double min_dis = core::MAX_D;
            for (size_t j = 0; j < Tree_nodes_fragile_V.size(); ++j) {
                // Tree_nodes_fragile_V[j] 预期是 3x1，分别存放 x,y,z
                const double x = Tree_nodes_fragile_V[j](0, 0);
                const double y = Tree_nodes_fragile_V[j](1, 0);
                const double z = Tree_nodes_fragile_V[j](2, 0);
                if (x >= boundary_left && x <= boundary_right &&
                    y >= boundary_bottom && y <= boundary_top) {
                    const double dz = std::abs(z - all_cut_layers[t][0].z);
                    if (dz < min_dis) {
                        min_dis = dz;
                    }
                }
            }

            if (min_dis <= threshold) {
                sum_value_of_fragile += 1.0;
            }
        }

        all_calculated_value.value_of_fragile = sum_value_of_fragile;
    }

    double calculate_projected_area(
        const Slicer& all_slicer,
        const std::vector<std::vector<vasco::core::Tri3>>& all_furcation_of_blocks,
        const std::vector<std::vector<cv::Point3d>>& all_cut_layers) noexcept
    {
        using Traits = CGAL::Convex_hull_traits_adapter_2<Kernel, CGAL::Pointer_property_map<Point_2>::type>;

        if (all_cut_layers.empty()) {
            return 0.0;
        }

        double avg_value = 0.0;

        for (size_t i = 0; i < all_cut_layers.size(); ++i) {
            // 界面面积
            std::vector<Point_2> points_in_interface;
            points_in_interface.reserve(all_cut_layers[i].size());
            for (const auto& p : all_cut_layers[i]) {
                points_in_interface.emplace_back(p.x, p.y);
            }
            double area_of_inter_face = 0.0;
            for (size_t j = 0; j < points_in_interface.size(); ++j) {
                const auto& a = points_in_interface[j];
                const auto& b = points_in_interface[(j + 1) % points_in_interface.size()];
                area_of_inter_face += (a.x() * b.y()) - (b.x() * a.y());
            }
            area_of_inter_face = std::abs(area_of_inter_face * 0.5);
            if (area_of_inter_face <= 0.0) {
                continue;
            }

            // 候选三角形投影点与凸包
            std::vector<Point_2> points_in_furcation_of_blocks;
            for (const auto& tri : all_furcation_of_blocks[i]) {
                for (int k = 0; k < 3; ++k) {
                    const auto vidx = tri[k];
                    points_in_furcation_of_blocks.emplace_back(
                        all_slicer.positions[vidx][0], all_slicer.positions[vidx][1]);
                }
            }

            std::vector<std::size_t> indices(points_in_furcation_of_blocks.size()), out;
            std::iota(indices.begin(), indices.end(), 0);
            CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
                Traits(CGAL::make_property_map(points_in_furcation_of_blocks)));

            std::vector<Point_2> points_in_projected_boundary;
            points_in_projected_boundary.reserve(out.size());
            for (std::size_t t : out) {
                points_in_projected_boundary.push_back(points_in_furcation_of_blocks[t]);
            }

            double area_of_projected_boundary = 0.0;
            for (size_t j = 0; j < points_in_projected_boundary.size(); ++j) {
                const auto& a = points_in_projected_boundary[j];
                const auto& b = points_in_projected_boundary[(j + 1) % points_in_projected_boundary.size()];
                area_of_projected_boundary += (a.x() * b.y()) - (b.x() * a.y());
            }
            area_of_projected_boundary = std::abs(area_of_projected_boundary * 0.5);

            // 多边形交叠面积的容斥计算（与原实现一致）
            PolyIntersec polyint;
            const int n1 = static_cast<int>(points_in_projected_boundary.size());
            const int n2 = static_cast<int>(points_in_interface.size());
            MyPoint ps1[10000];
            MyPoint ps2[10000];
            for (int j = 0; j < n1; ++j) {
                ps1[j].x = points_in_projected_boundary[j].x();
                ps1[j].y = points_in_projected_boundary[j].y();
            }
            for (int j = 0; j < n2; ++j) {
                ps2[j].x = points_in_interface[j].x();
                ps2[j].y = points_in_interface[j].y();
            }
            double inter = polyint.intersectArea(ps1, n1, ps2, n2);
            // 容斥：两者面积和减去交集
            double union_area = std::abs(polyint.area(ps1, n1)) + std::abs(polyint.area(ps2, n2)) - inter;

            // 使用原始的比例定义：凸包投影面积 / 切割界面面积
            avg_value += (area_of_projected_boundary / area_of_inter_face);
        }

        avg_value /= static_cast<double>(all_cut_layers.size());
        return avg_value;
    }

    std::vector<std::vector<int>> FindAllCutLayers(
        const Layer_Graph& layer_graph,
        const std::vector<std::vector<int>>& final_pathes,
        std::vector<std::vector<int>>& all_cut_layers_dependency_layer,
        bool& jud_admit) noexcept
    {
        std::vector<std::vector<int>> all_id_of_cut_layers;
        all_id_of_cut_layers.resize(final_pathes.size());
        all_cut_layers_dependency_layer.resize(final_pathes.size());

        for (size_t i = 0; i < final_pathes.size(); i++) {
            for (size_t j = 0; j < final_pathes[i].size(); j++) {
                bool jud_cut_layer = true;
                bool jud_cut_layer_fork = false;

				int layerID = final_pathes[i][j];

                if (layerID > 0) {
                    for (size_t k = 0; k < layer_graph.G_2[layerID].size(); k++) {
                        bool jud_cut_layer_2 = true;
                        for (size_t m = 0; m < final_pathes[i].size(); m++) {
                            if (j != m) {
								int fromID = layer_graph.edges[layer_graph.G_2[layerID][k]].GetFrom();
                                if (fromID == final_pathes[i][m]) {
                                    jud_cut_layer_2 = false;
                                    break;
                                }
                            }
                        }
                        if (!jud_cut_layer_2) {
                            jud_cut_layer = false;
                        }
                        if (jud_cut_layer_2) {
                            jud_cut_layer_fork = true;
                        }
                    }
                }
                

                // 对于中间存在分叉的情况，先剔除这些结果（保持与旧逻辑一致）
                if (!jud_cut_layer && jud_cut_layer_fork) {
                    jud_admit = false;
                    return all_id_of_cut_layers;
                }

                if (jud_cut_layer || layer_graph.G_2[layerID].empty()) {
                    // 删除孤立 layer
                    if (layer_graph.G_2[layerID].empty() && layer_graph.G[layerID].empty()) {
                        continue;
                    }
                    all_id_of_cut_layers[i].push_back(static_cast<int>(j));
                    all_cut_layers_dependency_layer[i].push_back(static_cast<int>(layer_graph.G_2[layerID].size()));
                }
            }
        }
        return all_id_of_cut_layers;
    }

} // namespace vasco


//void HybridManufacturing::SortEdges(vector<vector<area_S>> ori_all_the_area_S)
//{
//	for (int i = 0; i < ori_all_the_area_S.size(); i++) {
//		for (int j = 0; j < ori_all_the_area_S[i].size(); j++) {
//			for (int k = j + 1; k < ori_all_the_area_S[i].size(); k++) {
//				if (ori_all_the_area_S[i][j].oriId > ori_all_the_area_S[i][k].oriId)
//					swap(ori_all_the_area_S[i][j], ori_all_the_area_S[i][k]);
//			}
//		}
//	}
//}
//double HybridManufacturing::calculate_manufacturing_require(Slicer_2 slicer, int index, double min_z, double max_z)
//{
//    double average_z = (slicer.positions[slicer.triangles[index][0]][2] + slicer.positions[slicer.triangles[index][1]][2] + slicer.positions[slicer.triangles[index][2]][2]) / 3;
//    double value_of_manufacturing_require = (max_z - average_z) / (max_z - min_z);
//    return value_of_manufacturing_require;
//}


//void HybridManufacturing::Anticlockwise(vector<vector<int>>& real_cutting_plane_triangles, Slicer_2 all_slicer)
//{
//
//    for (int t = 0; t < real_cutting_plane_triangles.size(); t++) {
//        double max_x = -MAX_D;
//        int index_point;
//        for (int i = 0; i < real_cutting_plane_triangles[t].size(); i++)
//        {
//            if (all_slicer.positions[real_cutting_plane_triangles[t][i]][0] > max_x) {
//                max_x = all_slicer.positions[real_cutting_plane_triangles[t][i]][0];
//                index_point = i;
//            }
//        }
//        double d = (all_slicer.positions[real_cutting_plane_triangles[t][(index_point) % real_cutting_plane_triangles[t].size()]][0] - all_slicer.positions[real_cutting_plane_triangles[t][(index_point - 1 + real_cutting_plane_triangles[t].size()) % real_cutting_plane_triangles[t].size()]][0])
//            * (all_slicer.positions[real_cutting_plane_triangles[t][(index_point) % real_cutting_plane_triangles[t].size()]][1] - all_slicer.positions[real_cutting_plane_triangles[t][(index_point + 1) % real_cutting_plane_triangles[t].size()]][1])
//            - (all_slicer.positions[real_cutting_plane_triangles[t][(index_point) % real_cutting_plane_triangles[t].size()]][1] - all_slicer.positions[real_cutting_plane_triangles[t][(index_point - 1 + real_cutting_plane_triangles[t].size()) % real_cutting_plane_triangles[t].size()]][1])
//            * (all_slicer.positions[real_cutting_plane_triangles[t][(index_point) % real_cutting_plane_triangles[t].size()]][0] - all_slicer.positions[real_cutting_plane_triangles[t][(index_point + 1) % real_cutting_plane_triangles[t].size()]][0]);
//        if (d > 0)
//            reverse(real_cutting_plane_triangles[t].begin(), real_cutting_plane_triangles[t].end());
//    }
//
//
//    /*for (int t = 0; t < real_cutting_plane_triangles.size(); t++) {
//        double d = 0;
//        for (int i = 0; i < real_cutting_plane_triangles[t].size(); i++)
//        {
//            d += (all_slicer.positions[real_cutting_plane_triangles[t][(i + 1)% real_cutting_plane_triangles[t].size()]][0] - all_slicer.positions[real_cutting_plane_triangles[t][i]][0])
//                * (all_slicer.positions[real_cutting_plane_triangles[t][(i + 1)% real_cutting_plane_triangles[t].size()]][1] - all_slicer.positions[real_cutting_plane_triangles[t][(i + 2)% real_cutting_plane_triangles[t].size()]][1])
//                - (all_slicer.positions[real_cutting_plane_triangles[t][(i + 1)% real_cutting_plane_triangles[t].size()]][1] - all_slicer.positions[real_cutting_plane_triangles[t][i]][1])
//                * (all_slicer.positions[real_cutting_plane_triangles[t][(i + 1)% real_cutting_plane_triangles[t].size()]][0] - all_slicer.positions[real_cutting_plane_triangles[t][(i + 2)% real_cutting_plane_triangles[t].size()]][0]);
//        }
//        if (d > 0)
//            reverse(real_cutting_plane_triangles[t].begin(), real_cutting_plane_triangles[t].end());
//    }*/
//
//    /*for (int t = 0; t < real_cutting_plane_triangles.size(); t++) {
//        double d = 0;
//        for (int i = 0; i < real_cutting_plane_triangles[t].size() - 1; i++)
//        {
//            cv::Point2d temp_p1(all_slicer.positions[real_cutting_plane_triangles[t][i]][0], all_slicer.positions[real_cutting_plane_triangles[t][i]][1]);
//            cv::Point2d temp_p2(all_slicer.positions[real_cutting_plane_triangles[t][i+1]][0], all_slicer.positions[real_cutting_plane_triangles[t][i+1]][1]);
//            d += -0.5 * (all_slicer.positions[real_cutting_plane_triangles[t][i + 1]][1] + all_slicer.positions[real_cutting_plane_triangles[t][i]][1]) * (all_slicer.positions[real_cutting_plane_triangles[t][i + 1]][0] - all_slicer.positions[real_cutting_plane_triangles[t][i]][0]);
//        }
//        if (d < 0)
//            reverse(real_cutting_plane_triangles[t].begin(), real_cutting_plane_triangles[t].end());
//    }*/
//}

//void HybridManufacturing::calculate_fragile_value(all_value& all_calculated_value, vector<vector<cv::Point3d>> all_cut_layers, vector<Eigen::MatrixXd> Tree_nodes_fragile_V)
//{
//    double threshold = 6.0;
//    double sum_value_of_fragile = 0;
//    for (int t = 0; t < all_cut_layers.size(); t++) {
//        double boundary_bottom = 999999, boundary_left = 999999, boundary_top = -999999, boundary_right = -999999;
//        for (int i = 0; i < all_cut_layers[t].size(); i++) {
//            boundary_top = std::max(boundary_top, all_cut_layers[t][i].y);
//            boundary_bottom = std::min(boundary_bottom, all_cut_layers[t][i].y);
//            boundary_right = std::max(boundary_right, all_cut_layers[t][i].x);
//            boundary_left = std::min(boundary_left, all_cut_layers[t][i].x);
//        }
//        double min_dis = MAX_D;
//        for (int j = 0; j < Tree_nodes_fragile_V.size(); j++) {
//            if (Tree_nodes_fragile_V[j](0, 0) >= boundary_left && Tree_nodes_fragile_V[j](0, 0) <= boundary_right && Tree_nodes_fragile_V[j](1, 0) >= boundary_bottom && Tree_nodes_fragile_V[j](1, 0) <= boundary_top) {
//                if (abs(Tree_nodes_fragile_V[j](2, 0) - all_cut_layers[t][0].z) < min_dis)
//                    min_dis = abs(Tree_nodes_fragile_V[j](2, 0) - all_cut_layers[t][0].z);
//            }
//        }
//
//        if (min_dis <= threshold)
//            sum_value_of_fragile += 1;
//    }
//    all_calculated_value.value_of_fragile = sum_value_of_fragile;
//}

//double HybridManufacturing::calculate_projected_area(Slicer_2 all_slicer, vector<vector<TRiangle>> all_furcation_of_blocks, vector<vector<cv::Point3d>> all_cut_layers)
//{
//    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//    typedef K::Point_2 Point_2;
//    typedef CGAL::Convex_hull_traits_adapter_2<K,
//        CGAL::Pointer_property_map<Point_2>::type > Convex_hull_traits_2;
//    vector<vector<cv::Point3d>> vis_points;
//    vis_points.resize(all_cut_layers.size());
//    double avg_value = 0;
//
//    for (int i = 0; i < all_cut_layers.size(); i++) {
//        vector<Point_2> points_in_interface;
//        points_in_interface.clear();
//        double area_of_inter_face = 0;
//        double area_of_projected_boundary = 0;
//        vector<Point_2> points_in_furcation_of_blocks;
//        points_in_furcation_of_blocks.clear();
//        vector<Point_2> points_in_projected_boundary;
//        points_in_projected_boundary.clear();
//
//        for (int j = 0; j < all_cut_layers[i].size(); j++) {
//            Point_2 temp_point(all_cut_layers[i][j].x, all_cut_layers[i][j].y);
//            points_in_interface.push_back(temp_point);
//        }
//        for (int j = 0; j < points_in_interface.size(); j++) {
//            area_of_inter_face += (points_in_interface[j].x() * points_in_interface[(j + 1) % points_in_interface.size()].y()) - (points_in_interface[(j + 1) % points_in_interface.size()].x() * points_in_interface[j].y());
//        }
//        area_of_inter_face /= 2;
//        area_of_inter_face = abs(area_of_inter_face);
//
//        for (int j = 0; j < all_furcation_of_blocks[i].size(); j++) {
//            for (int k = 0; k < 3; k++) {
//                Point_2 temp_point(all_slicer.positions[all_furcation_of_blocks[i][j][k]][0], all_slicer.positions[all_furcation_of_blocks[i][j][k]][1]);
//                points_in_furcation_of_blocks.push_back(temp_point);
//                //cv::Point3d temp_3d_point(temp_point.x(), temp_point.y(), 0);
//                //vis_points[i].push_back(temp_3d_point);
//            }
//        }
//        std::vector<std::size_t> indices(points_in_furcation_of_blocks.size()), out;
//        std::iota(indices.begin(), indices.end(), 0);
//        CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
//            Convex_hull_traits_2(CGAL::make_property_map(points_in_furcation_of_blocks)));
//        for (std::size_t t : out) {
//            //std::cout << "points[" << t << "] = " << points_in_furcation_of_blocks[t] << std::endl;
//            points_in_projected_boundary.push_back(points_in_furcation_of_blocks[t]);
//            //cv::Point3d temp_3d_point(points_in_furcation_of_blocks[t].x(), points_in_furcation_of_blocks[t].y(), 0);
//            //vis_points[i].push_back(temp_3d_point);
//        }
//
//        for (int j = 0; j < points_in_projected_boundary.size(); j++) {
//            area_of_projected_boundary += (points_in_projected_boundary[j].x() * points_in_projected_boundary[(j + 1) % points_in_projected_boundary.size()].y()) - (points_in_projected_boundary[(j + 1) % points_in_projected_boundary.size()].x() * points_in_projected_boundary[j].y());
//        }
//        area_of_projected_boundary /= 2;
//        area_of_projected_boundary = abs(area_of_projected_boundary);
//
//
//        int n1 = points_in_projected_boundary.size(), n2 = points_in_interface.size();
//        PolyIntersec polyint;
//        MyPoint ps1[10000];
//        MyPoint ps2[10000];
//        for (int j = 0; j < n1; j++) {
//            ps1[j].x = points_in_projected_boundary[j].x();
//            ps1[j].y = points_in_projected_boundary[j].y();
//        }
//        for (int j = 0; j < n2; j++) {
//            ps2[j].x = points_in_interface[j].x();
//            ps2[j].y = points_in_interface[j].y();
//        }
//        double ans = polyint.intersectArea(ps1, n1, ps2, n2);
//        ans = fabs(polyint.area(ps1, n1)) + fabs(polyint.area(ps2, n2)) - ans;//容斥
//        //cout << ans<<endl;
//
//        avg_value += ((area_of_projected_boundary) / area_of_inter_face);
//    }
//    avg_value /= all_cut_layers.size();
//
//    //Visual vis;
//    //vis.generateModelForRendering_9(vis_points, ".\\vis\\projected");
//    return avg_value;
//}