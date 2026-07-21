#pragma once
#ifndef SURFACE_MESH_SLICE_DATA_H_
#define SURFACE_MESH_SLICE_DATA_H_

#include <vector>
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Point_2 = Kernel::Point_2;
using SurfaceMesh = CGAL::Surface_mesh<Point_3>;

using Polyline_type = std::vector<Point_3>;
using Polylines = std::vector<Polyline_type>;

struct SurfaceMeshSliceSegment {
	Point_3 start;
	Point_3 end;
	int face_id = -1;
};

struct SurfaceMeshSliceData {
	Polylines contour_points;
	// Local contour index of the direct parent in this slice; -1 means root.
	std::vector<int> contour_parent_ids;
	std::vector<std::vector<SurfaceMesh::Face_index>> contour_face_ids;
	double layer_z = 0.0;
	std::vector<Eigen::Vector3d> face_normals;
};

#endif
