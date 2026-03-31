#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h> // bounded_side_2
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>

// 统一使用 CGAL Kernel 类型
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

namespace vasco {


/*
 * 检测点 pt 与多边形范围关系：
 * 返回 true: 位于有界内部或边界上
 * 返回 false: 位于无界外部
 * 说明：
 *  - pgn_begin / pgn_end 为 ForwardIterator，解引用后应为 Point_2。
 *  - traits 传入 Kernel() 或已有的 Kernel 实例。
 */
template <class ForwardIterator>
inline bool check_inside_2(const Point_2& pt,
                           ForwardIterator pgn_begin,
                           ForwardIterator pgn_end,
                           const Kernel& traits = Kernel())
{
    switch (CGAL::bounded_side_2(pgn_begin, pgn_end, pt, traits)) {
    case CGAL::ON_BOUNDED_SIDE:
    case CGAL::ON_BOUNDARY:
        return true;
    case CGAL::ON_UNBOUNDED_SIDE:
        return false;
    default:
        return false;
    }
}

} // namespace vasco