#include "vasco/core/ContourContainment.h"

#include <CGAL/Polygon_2.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace vasco::contour_containment
{
	namespace
	{
		using ContainmentKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
		using ContainmentPoint2 = ContainmentKernel::Point_2;
		using ContainmentPolygon2 = CGAL::Polygon_2<ContainmentKernel>;

		struct PreparedContour {
			ContainmentPolygon2 polygon;
			double absolute_area = 0.0;
			double min_x = 0.0;
			double min_y = 0.0;
			double max_x = 0.0;
			double max_y = 0.0;
			bool valid = false;
		};

		// Computes planar squared distance without an unnecessary square root.
		double SquaredDistance(const ContainmentPoint2& a, const ContainmentPoint2& b)
		{
			const double dx = CGAL::to_double(a.x() - b.x());
			const double dy = CGAL::to_double(a.y() - b.y());
			return dx * dx + dy * dy;
		}

		// Cleans one 3D slice contour and prepares its 2D polygon, area, and AABB.
		// Invalid, degenerate, or self-intersecting contours remain marked invalid.
		PreparedContour PrepareContour(const Polyline_type& contour, double tolerance)
		{
			PreparedContour prepared;
			std::vector<ContainmentPoint2> points;
			points.reserve(contour.size());
			const double tolerance_squared = tolerance * tolerance;

			for (const Point_3& point : contour) {
				const double x = CGAL::to_double(point.x());
				const double y = CGAL::to_double(point.y());
				if (!std::isfinite(x) || !std::isfinite(y)) {
					continue;
				}

				const ContainmentPoint2 point_2(x, y);
				if (!points.empty() && SquaredDistance(points.back(), point_2) <= tolerance_squared) {
					continue;
				}
				points.push_back(point_2);
			}

			if (points.size() >= 2
				&& SquaredDistance(points.front(), points.back()) <= tolerance_squared) {
				points.pop_back();
			}
			if (points.size() < 3) {
				return prepared;
			}

			prepared.polygon = ContainmentPolygon2(points.begin(), points.end());
			if (!prepared.polygon.is_simple()) {
				return prepared;
			}

			prepared.absolute_area = std::abs(CGAL::to_double(prepared.polygon.area()));
			if (!std::isfinite(prepared.absolute_area) || prepared.absolute_area <= tolerance_squared) {
				return prepared;
			}

			prepared.min_x = prepared.max_x = CGAL::to_double(points.front().x());
			prepared.min_y = prepared.max_y = CGAL::to_double(points.front().y());
			for (const ContainmentPoint2& point : points) {
				const double x = CGAL::to_double(point.x());
				const double y = CGAL::to_double(point.y());
				prepared.min_x = std::min(prepared.min_x, x);
				prepared.min_y = std::min(prepared.min_y, y);
				prepared.max_x = std::max(prepared.max_x, x);
				prepared.max_y = std::max(prepared.max_y, y);
			}
			prepared.valid = true;
			return prepared;
		}

		// Performs a cheap AABB containment test before the polygon-side predicate.
		bool BoundingBoxContains(
			const PreparedContour& parent,
			const PreparedContour& child,
			double tolerance)
		{
			return parent.min_x <= child.min_x + tolerance
				&& parent.min_y <= child.min_y + tolerance
				&& parent.max_x + tolerance >= child.max_x
				&& parent.max_y + tolerance >= child.max_y;
		}
	}

	// Selects the smallest valid enclosing contour as each node's direct parent.
	// Area and AABB filters avoid most expensive point-in-polygon evaluations.
	std::vector<int> BuildContourParentIds(
		const Polylines& contours,
		double coordinate_tolerance)
	{
		const double tolerance = std::max(0.0, coordinate_tolerance);
		std::vector<PreparedContour> prepared_contours;
		prepared_contours.reserve(contours.size());
		for (const Polyline_type& contour : contours) {
			prepared_contours.push_back(PrepareContour(contour, tolerance));
		}

		std::vector<int> parent_ids(contours.size(), -1);
		for (int child_id = 0; child_id < static_cast<int>(prepared_contours.size()); ++child_id) {
			const PreparedContour& child = prepared_contours[child_id];
			if (!child.valid || child.polygon.size() == 0) {
				continue;
			}

			double smallest_parent_area = std::numeric_limits<double>::max();
			for (int candidate_id = 0; candidate_id < static_cast<int>(prepared_contours.size()); ++candidate_id) {
				if (candidate_id == child_id) {
					continue;
				}

				const PreparedContour& candidate = prepared_contours[candidate_id];
				if (!candidate.valid
					|| candidate.absolute_area <= child.absolute_area
					|| candidate.absolute_area >= smallest_parent_area
					|| !BoundingBoxContains(candidate, child, tolerance)) {
					continue;
				}

				const auto side = candidate.polygon.bounded_side(child.polygon[0]);
				if (side != CGAL::ON_BOUNDED_SIDE) {
					continue;
				}

				parent_ids[child_id] = candidate_id;
				smallest_parent_area = candidate.absolute_area;
			}
		}
		return parent_ids;
	}

	// Reverses the parent-only representation into direct-child adjacency lists.
	std::vector<std::vector<int>> BuildDirectChildren(
		const std::vector<int>& parent_ids)
	{
		std::vector<std::vector<int>> children(parent_ids.size());
		for (int child_id = 0; child_id < static_cast<int>(parent_ids.size()); ++child_id) {
			const int parent_id = parent_ids[child_id];
			if (parent_id >= 0
				&& parent_id < static_cast<int>(parent_ids.size())
				&& parent_id != child_id) {
				children[parent_id].push_back(child_id);
			}
		}
		return children;
	}

	// Walks a node's parent chain to determine nesting depth and detect cycles.
	int ComputeContourDepth(
		const std::vector<int>& parent_ids,
		int contour_id)
	{
		if (contour_id < 0 || contour_id >= static_cast<int>(parent_ids.size())) {
			return -1;
		}

		std::vector<bool> visited(parent_ids.size(), false);
		int depth = 0;
		int current_id = contour_id;
		while (parent_ids[current_id] != -1) {
			if (visited[current_id]) {
				return -1;
			}
			visited[current_id] = true;

			current_id = parent_ids[current_id];
			if (current_id < 0 || current_id >= static_cast<int>(parent_ids.size())) {
				return -1;
			}
			++depth;
		}
		return depth;
	}

	// Converts the parent-only tree into material regions using even-odd depth.
	std::vector<ContourComponent> BuildMaterialComponents(
		const std::vector<int>& parent_ids)
	{
		const std::vector<std::vector<int>> children = BuildDirectChildren(parent_ids);
		std::vector<ContourComponent> components;
		components.reserve(parent_ids.size());
		for (int contour_id = 0; contour_id < static_cast<int>(parent_ids.size()); ++contour_id) {
			const int depth = ComputeContourDepth(parent_ids, contour_id);
			if (depth < 0 || depth % 2 != 0) {
				continue;
			}

			ContourComponent component;
			component.outer_contour_id = contour_id;
			for (int child_id : children[contour_id]) {
				const int child_depth = ComputeContourDepth(parent_ids, child_id);
				if (child_depth == depth + 1 && child_depth % 2 != 0) {
					component.hole_contour_ids.push_back(child_id);
				}
			}
			components.push_back(std::move(component));
		}
		return components;
	}
}
