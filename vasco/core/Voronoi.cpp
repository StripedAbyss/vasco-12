#include "Voronoi.h"
#include <Eigen/Geometry>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/centroid.h>
#include <iostream>

namespace vasco
{
	void BuildVoronoiCells(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		double thresholdZ,
		std::vector<VoronoiCell>& outCells,
		std::vector<Eigen::Vector3d>& outBottomVertices,
		bool visualize,
		const std::string& fileNameNoExt)
	{
		outCells.clear();
		outBottomVertices.clear();

		int start_time = clock();
		std::vector<std::vector<Eigen::Vector3d>> all_lines;
		double min_z = 1e100; //所有顶点z值的最小值
		for (int i = 0; i < V.rows(); ++i) {
			min_z = std::min(min_z, V(i, 2));
		}

		auto collect_adjacent_faces = [&](int vertex_index) {
			std::vector<int> adjacent_faces;
			for (int j = 0; j < F.rows(); ++j) {
				for (int k = 0; k < 3; ++k) {
					if (F(j, k) == vertex_index) {
						adjacent_faces.push_back(j);
					}
				}
			}
			return adjacent_faces;
			};

		for (int i = 0; i < V.rows(); ++i)
		{
			if (V(i, 2) - min_z <= thresholdZ)
			{
				outBottomVertices.emplace_back(V(i, 0), V(i, 1), V(i, 2));
				VoronoiCell new_cell;
				new_cell.is_available = false;
				outCells.push_back(new_cell);
				continue;
			}

			VoronoiCell new_cell;
			std::vector<Eigen::MatrixXd> all_boundary_V;
			std::vector<int> index_of_adjacent_F = collect_adjacent_faces(i);
			std::vector<int> index_of_adjacent_V;

			int index_begin = -1;
			int index_end = -1;
			int index_first_begin = -1;
			for (int j = 0; j < 3; ++j)
			{
				if (F(index_of_adjacent_F[0], j) == i)
				{
					index_begin = F(index_of_adjacent_F[0], (j + 1) % 3);
					index_end = F(index_of_adjacent_F[0], (j + 2) % 3);
					index_first_begin = index_begin;
					break;
				}
			}
			index_of_adjacent_V.push_back(index_begin);
			int cont_num = 1;
			while (index_end != index_first_begin)
			{
				for (int j = cont_num; j < static_cast<int>(index_of_adjacent_F.size()); ++j)
				{
					for (int k = 0; k < 3; ++k)
					{
						if (F(index_of_adjacent_F[j], k) != i)
						{
							continue;
						}
						if (F(index_of_adjacent_F[j], (k + 1) % 3) == index_end)
						{
							index_begin = F(index_of_adjacent_F[j], (k + 1) % 3);
							index_end = F(index_of_adjacent_F[j], (k + 2) % 3);
							index_of_adjacent_V.push_back(index_begin);
							std::swap(index_of_adjacent_F[j], index_of_adjacent_F[cont_num]);
							++cont_num;
						}
						else if (F(index_of_adjacent_F[j], (k + 2) % 3) == index_end)
						{
							index_begin = F(index_of_adjacent_F[j], (k + 2) % 3);
							index_end = F(index_of_adjacent_F[j], (k + 1) % 3);
							index_of_adjacent_V.push_back(index_begin);
							std::swap(index_of_adjacent_F[j], index_of_adjacent_F[cont_num]);
							++cont_num;
						}
					}
				}
			}

			for (int j = 0; j < (int)index_of_adjacent_F.size(); ++j)
			{
				double na = (V(F(index_of_adjacent_F[j], 1), 1) - V(F(index_of_adjacent_F[j], 0), 1)) *
					(V(F(index_of_adjacent_F[j], 2), 2) - V(F(index_of_adjacent_F[j], 0), 2)) -
					(V(F(index_of_adjacent_F[j], 1), 2) - V(F(index_of_adjacent_F[j], 0), 2)) *
					(V(F(index_of_adjacent_F[j], 2), 1) - V(F(index_of_adjacent_F[j], 0), 1));
				double nb = (V(F(index_of_adjacent_F[j], 1), 2) - V(F(index_of_adjacent_F[j], 0), 2)) *
					(V(F(index_of_adjacent_F[j], 2), 0) - V(F(index_of_adjacent_F[j], 0), 0)) -
					(V(F(index_of_adjacent_F[j], 1), 0) - V(F(index_of_adjacent_F[j], 0), 0)) *
					(V(F(index_of_adjacent_F[j], 2), 2) - V(F(index_of_adjacent_F[j], 0), 2));
				double nc = (V(F(index_of_adjacent_F[j], 1), 0) - V(F(index_of_adjacent_F[j], 0), 0)) *
					(V(F(index_of_adjacent_F[j], 2), 1) - V(F(index_of_adjacent_F[j], 0), 1)) -
					(V(F(index_of_adjacent_F[j], 1), 1) - V(F(index_of_adjacent_F[j], 0), 1)) *
					(V(F(index_of_adjacent_F[j], 2), 0) - V(F(index_of_adjacent_F[j], 0), 0));
				Eigen::Vector3d vn(na, nb, nc);

				Eigen::Vector3d vectorBefore(0, 0, 1);
				Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vn).toRotationMatrix();

				std::vector<Eigen::MatrixXd> current_V;
				current_V.reserve(3);
				Eigen::MatrixXd temp_V(3, 1);

				temp_V << V(i, 0), V(i, 1), V(i, 2);
				current_V.push_back(temp_V);
				temp_V << V(index_of_adjacent_V[j], 0), V(index_of_adjacent_V[j], 1), V(index_of_adjacent_V[j], 2);
				current_V.push_back(temp_V);
				temp_V << V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 0),
					V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 1),
					V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 2);
				current_V.push_back(temp_V);

				for (auto& mv : current_V)
					mv = rotMatrix.inverse() * mv;

				// ʹ������
				Eigen::MatrixXd new_V(3, 1);
				new_V(0, 0) = (current_V[0](0, 0) + current_V[1](0, 0) + current_V[2](0, 0)) / 3.0;
				new_V(1, 0) = (current_V[0](1, 0) + current_V[1](1, 0) + current_V[2](1, 0)) / 3.0;
				new_V(2, 0) = current_V[0](2, 0);

				rotMatrix = Eigen::Quaterniond::FromTwoVectors(vn, vectorBefore).toRotationMatrix();
				new_V = rotMatrix.inverse() * new_V;

				all_boundary_V.push_back(new_V);

				if (all_boundary_V.size() == 2)
				{
					Eigen::Vector3d v1(V(i, 0), V(i, 1), V(i, 2));
					Eigen::Vector3d v2(all_boundary_V[0](0, 0), all_boundary_V[0](1, 0), all_boundary_V[0](2, 0));
					Eigen::Vector3d v3(all_boundary_V[1](0, 0), all_boundary_V[1](1, 0), all_boundary_V[1](2, 0));
					double na2 = (v2.y() - v1.y()) * (v3.z() - v1.z()) - (v2.z() - v1.z()) * (v3.y() - v1.y());
					double nb2 = (v2.z() - v1.z()) * (v3.x() - v1.x()) - (v2.x() - v1.x()) * (v3.z() - v1.z());
					double nc2 = (v2.x() - v1.x()) * (v3.y() - v1.y()) - (v2.y() - v1.y()) * (v3.x() - v1.x());
					Eigen::Vector3d vn2(na2, nb2, nc2);
					if (vn.dot(vn2) < 0)
						std::swap(all_boundary_V[0], all_boundary_V[1]);
				}
			}


			new_cell.is_available = true;
			new_cell.site = i;
			new_cell.adjacent_cells = index_of_adjacent_V;

			std::vector<Eigen::Vector3d> temp_lines;
			all_lines.push_back(temp_lines);
			for (const auto& bV : all_boundary_V)
			{
				Eigen::Vector3d temp_vec(bV(0, 0), bV(1, 0), bV(2, 0));
				all_lines.back().push_back(temp_vec);
				new_cell.all_points_in_polygon.emplace_back(temp_vec.x(), temp_vec.y(), temp_vec.z());
			}
			outCells.push_back(std::move(new_cell));
		}

		int end_time = clock();
		std::cout << "&&&time&&& Voronoi build: " << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;

		if (visualize)
		{
			Visual vis;
			vis.generateModelForRendering_8(all_lines, fileNameNoExt + "_voronoi.obj");
			std::cout << "visualized " << fileNameNoExt + "_voronoi.obj" << std::endl;
		}
	}

	void BuildVoronoiCells(const SurfaceMesh& mesh,
		double thresholdZ,
		std::vector<VoronoiCell>& outCells,
		std::vector<Eigen::Vector3d>& outBottomVertices,
		bool visualize,
		const std::string& fileNameNoExt)
	{
		outCells.clear();
		outBottomVertices.clear();

		int start_time = clock();
		std::vector<std::vector<Eigen::Vector3d>> all_lines;

		double min_z = 1e100;
		std::vector<SurfaceMesh::Vertex_index> vertices;
		vertices.reserve(mesh.number_of_vertices());
		int v_cnt = 0;
		for (auto v : mesh.vertices()) {
			vertices.push_back(v);
			const auto& p = mesh.point(v);
			if (v.id() != v_cnt) {
				std::cerr << "Error: vertex index does not match its id!" << std::endl;
				return;
			}
			min_z = std::min(min_z, p.z());
			v_cnt++;
		}

		for (int i = 0; i < static_cast<int>(vertices.size()); ++i)
		{
			const auto vi = vertices[i];
			const auto& vi_point = mesh.point(vi);
			if (vi_point.z() - min_z <= thresholdZ)
			{
				outBottomVertices.emplace_back(vi_point.x(), vi_point.y(), vi_point.z());
				VoronoiCell new_cell;
				new_cell.is_available = false;
				outCells.push_back(new_cell);
				continue;
			}

			std::vector<SurfaceMesh::Halfedge_index> around_halfedges;
			const auto start_halfedge = mesh.halfedge(vi);
			if (start_halfedge == SurfaceMesh::null_halfedge()) { //特判孤立点
				VoronoiCell new_cell;
				new_cell.is_available = false;
				outCells.push_back(new_cell);
				std::cout << "Warning: vertex " << i << " is isolated, no adjacent faces!" << std::endl;
				continue;
			}
			for (auto h : mesh.halfedges_around_target(start_halfedge)) {
				if (mesh.is_border(h)) {
					std::cout << "Warning: vertex " << i << " is adjacent to border halfedge, skipping this halfedge!" << std::endl;
					continue;
				}
				around_halfedges.push_back(h);
			}

			if (around_halfedges.empty()) { //特判相邻面为空的点
				VoronoiCell new_cell;
				new_cell.is_available = false;
				outCells.push_back(new_cell);
				std::cout << "Warning: vertex " << i << " has no adjacent faces!" << std::endl;
				continue;
			}

			std::vector<int> index_of_adjacent_V;
			index_of_adjacent_V.reserve(around_halfedges.size());
			for (const auto h : around_halfedges) {
				index_of_adjacent_V.push_back(static_cast<int>(mesh.source(h).idx()));
			}

			VoronoiCell new_cell;
			std::vector<Eigen::MatrixXd> all_boundary_V;
			all_boundary_V.reserve(around_halfedges.size());

			for (int j = 0; j < static_cast<int>(around_halfedges.size()); ++j)
			{
				auto h = around_halfedges[j];
				auto v_prev = mesh.source(h);
				auto v_next = mesh.target(mesh.next(h));
				const auto& p0 = mesh.point(v_prev);
				const auto& p1 = vi_point;
				const auto& p2 = mesh.point(v_next);

				const auto centroid = CGAL::centroid(p1, p0, p2);
				Eigen::Vector3d centroid_vec(centroid.x(), centroid.y(), centroid.z());

				Eigen::MatrixXd new_V(3, 1);
				new_V(0, 0) = centroid_vec.x();
				new_V(1, 0) = centroid_vec.y();
				new_V(2, 0) = centroid_vec.z();

				all_boundary_V.push_back(new_V);

			}

			new_cell.is_available = true;
			new_cell.site = static_cast<int>(vi.idx());
			new_cell.adjacent_cells = index_of_adjacent_V;

			std::vector<Eigen::Vector3d> temp_lines;
			all_lines.push_back(temp_lines);
			for (const auto& bV : all_boundary_V)
			{
				Eigen::Vector3d temp_vec(bV(0, 0), bV(1, 0), bV(2, 0));
				all_lines.back().push_back(temp_vec);
				new_cell.all_points_in_polygon.emplace_back(temp_vec.x(), temp_vec.y(), temp_vec.z());
			}
			outCells.push_back(std::move(new_cell));
		}

		int end_time = clock();
		std::cout << "&&&time&&& Voronoi build: " << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;

		if (visualize)
		{
			Visual vis;
			vis.generateModelForRendering_8(all_lines, fileNameNoExt + "_voronoi_zxx.obj");
			std::cout << "visualized " << fileNameNoExt + "_voronoi_zxx.obj" << std::endl;
		}
	}
}


//void HybridManufacturing::GetVoronoiCells()
//{
//	int start_time = clock();
//	double threshold_z = 2;  //�ײ��߶���ֵ
//	vector<vector<Eigen::Vector3d>> all_lines;
//	double min_z = MAX_D;
//	for (int i = 0; i < V.rows(); i++)
//		min_z = min(min_z, V(i, 2)); //�ҵ�������z�������Сֵ
//	for (int i = 0; i < V.rows(); i++) {
//		if (V(i, 2) - min_z > threshold_z) { //ö�����и�����ֵ�Ķ���i
//			VoronoiCell new_cell;	//ÿ��ö�ٵĶ��㴴���µ�Voronoi��Ԫ
//			//new_cell.site_normal.m_x = new_cell.site_normal.m_y = new_cell.site_normal.m_z = 0;
//			vector<Eigen::MatrixXd> all_boundary_V;	//�洢�õ�Ԫvoronoi������б߽綥��
//			vector<int> index_of_adjacent_F;
//			index_of_adjacent_F.clear();
//			vector<int> index_of_adjacent_V;
//			index_of_adjacent_V.clear();
//			for (int j = 0; j < F.rows(); j++) {
//				for (int k = 0; k < 3; k++) {
//					if (F(j, k) == i) {
//						index_of_adjacent_F.push_back(j); // �ҵ����а����ö���i����
//					}
//				}
//			}
//			int index_begin, index_end;
//			int index_first_begin;
//			for (int j = 0; j < 3; j++)
//				if (F(index_of_adjacent_F[0], j) == i) {
//					index_begin = F(index_of_adjacent_F[0], (j + 1) % 3); //�ҵ���һ����������������㣬��Ϊ�߽�����������յ�
//					index_end = F(index_of_adjacent_F[0], (j + 2) % 3);
//					index_first_begin = index_begin;
//					break;
//				}
//			index_of_adjacent_V.push_back(index_begin);	//���������ڽӶ����б�
//			int cont_num = 1;
//			while (index_end != index_first_begin) {
//				for (int j = cont_num; j < index_of_adjacent_F.size(); j++) {
//					for (int k = 0; k < 3; k++) {
//						if (F(index_of_adjacent_F[j], k) == i) {
//							if (F(index_of_adjacent_F[j], (k + 1) % 3) == index_end) {
//								index_begin = F(index_of_adjacent_F[j], (k + 1) % 3);	//�ҵ���һ���ڽ��棬������Ӧ�Ķ�������ڽӶ����б�
//								index_end = F(index_of_adjacent_F[j], (k + 2) % 3);
//								index_of_adjacent_V.push_back(index_begin);
//								swap(index_of_adjacent_F[j], index_of_adjacent_F[cont_num]);	//���ڽ��潻�������б���ǰ��
//								cont_num++;
//							}
//							else if (F(index_of_adjacent_F[j], (k + 2) % 3) == index_end) {
//								index_begin = F(index_of_adjacent_F[j], (k + 2) % 3);
//								index_end = F(index_of_adjacent_F[j], (k + 1) % 3);
//								index_of_adjacent_V.push_back(index_begin);
//								swap(index_of_adjacent_F[j], index_of_adjacent_F[cont_num]);
//								cont_num++;
//							}
//						}
//					}
//
//				}
//			}
//			for (int j = 0; j < index_of_adjacent_F.size(); j++) {
//				double na = (V(F(index_of_adjacent_F[j], 1), 1) - V(F(index_of_adjacent_F[j], 0), 1)) * (V(F(index_of_adjacent_F[j], 2), 2) - V(F(index_of_adjacent_F[j], 0), 2)) - (V(F(index_of_adjacent_F[j], 1), 2) - V(F(index_of_adjacent_F[j], 0), 2)) * (V(F(index_of_adjacent_F[j], 2), 1) - V(F(index_of_adjacent_F[j], 0), 1));
//				double nb = (V(F(index_of_adjacent_F[j], 1), 2) - V(F(index_of_adjacent_F[j], 0), 2)) * (V(F(index_of_adjacent_F[j], 2), 0) - V(F(index_of_adjacent_F[j], 0), 0)) - (V(F(index_of_adjacent_F[j], 1), 0) - V(F(index_of_adjacent_F[j], 0), 0)) * (V(F(index_of_adjacent_F[j], 2), 2) - V(F(index_of_adjacent_F[j], 0), 2));
//				double nc = (V(F(index_of_adjacent_F[j], 1), 0) - V(F(index_of_adjacent_F[j], 0), 0)) * (V(F(index_of_adjacent_F[j], 2), 1) - V(F(index_of_adjacent_F[j], 0), 1)) - (V(F(index_of_adjacent_F[j], 1), 1) - V(F(index_of_adjacent_F[j], 0), 1)) * (V(F(index_of_adjacent_F[j], 2), 0) - V(F(index_of_adjacent_F[j], 0), 0));
//				Eigen::Vector3d vn(na, nb, nc); //������������������ķ�����
//				//new_cell.site_normal.m_x += vn[0];
//				//new_cell.site_normal.m_y += vn[1];
//				//new_cell.site_normal.m_z += vn[2];
//				Eigen::Matrix3d rotMatrix;
//
//				Eigen::Vector3d vectorBefore(0, 0, 1);
//				rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vn).toRotationMatrix(); //�����z�ᵽ����������ת����
//				vector<Eigen::MatrixXd> current_V;	//�洢��ǰ�����������
//				Eigen::MatrixXd temp_V;
//				temp_V.resize(3, 1);
//				temp_V(0, 0) = V(i, 0);
//				temp_V(1, 0) = V(i, 1);
//				temp_V(2, 0) = V(i, 2);
//				current_V.push_back(temp_V);
//				temp_V(0, 0) = V(index_of_adjacent_V[j], 0);
//				temp_V(1, 0) = V(index_of_adjacent_V[j], 1);
//				temp_V(2, 0) = V(index_of_adjacent_V[j], 2);
//				current_V.push_back(temp_V);
//				temp_V(0, 0) = V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 0);
//				temp_V(1, 0) = V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 1);
//				temp_V(2, 0) = V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 2);
//				current_V.push_back(temp_V);
//				for (int k = 0; k < current_V.size(); k++) {
//					current_V[k] = rotMatrix.inverse() * current_V[k];	//�ѵ�ǰ��current_V������������ת���Է�����Ϊz�������ϵ��
//				}
//
//				double x1 = current_V[0](0, 0), x2 = current_V[1](0, 0);
//				double y1 = current_V[0](1, 0), y2 = current_V[1](1, 0);
//				double B1 = 1, A1 = (x1 - x2) / (y1 - y2); double C1 = -((y1 + y2) / 2 + A1 * (x1 + x2) / 2);	//���㴹ֱƽ���߲���
//				x1 = current_V[0](0, 0), x2 = current_V[2](0, 0);
//				y1 = current_V[0](1, 0), y2 = current_V[2](1, 0);
//				double B2 = 1, A2 = (x1 - x2) / (y1 - y2); double C2 = -((y1 + y2) / 2 + A2 * (x1 + x2) / 2);
//
//				Eigen::MatrixXd new_V;
//				new_V.resize(3, 1);
//				new_V(0, 0) = (C2 * B1 - C1 * B2) / (A1 * B2 - A2 * B1);	//���㴹ֱƽ���ߵĽ��㣬��Voronoi����
//				new_V(1, 0) = (C1 * A2 - C2 * A1) / (A1 * B2 - A2 * B1);
//				new_V(2, 0) = current_V[0](2, 0);
//
//				rotMatrix = Eigen::Quaterniond::FromTwoVectors(vn, vectorBefore).toRotationMatrix();
//				new_V = rotMatrix.inverse() * new_V;
//				all_boundary_V.push_back(new_V);	//������õ���Voronoi�������voronoi��߽綥���б�
//
//				if (all_boundary_V.size() == 2) {	//��֤�߽綥���˳��ʹ�ô�siteָ���һ���߽綥�㣬��ָ��ڶ����߽綥��ʱ��Ϊ��ʱ�뷽��
//					Eigen::Vector3d v1 = Eigen::Vector3d(V(i, 0), V(i, 1), V(i, 2));
//					Eigen::Vector3d v2 = Eigen::Vector3d(all_boundary_V[0](0, 0), all_boundary_V[0](1, 0), all_boundary_V[0](2, 0));
//					Eigen::Vector3d v3 = Eigen::Vector3d(all_boundary_V[1](0, 0), all_boundary_V[1](1, 0), all_boundary_V[1](2, 0));
//					//double ans = (v2.x() - v1.x()) * (v2.y() - v3.y()) - (v2.y() - v1.y()) * (v2.x() - v3.x());
//					//if (ans > 0)	//is clockwise
//					//	swap(v2, v3);
//					double na = (v2.y() - v1.y()) * (v3.z() - v1.z()) - (v2.z() - v1.z()) * (v3.y() - v1.y());
//					double nb = (v2.z() - v1.z()) * (v3.x() - v1.x()) - (v2.x() - v1.x()) * (v3.z() - v1.z());
//					double nc = (v2.x() - v1.x()) * (v3.y() - v1.y()) - (v2.y() - v1.y()) * (v3.x() - v1.x());
//					Eigen::Vector3d vn_2(na, nb, nc);
//					if (vn.dot(vn_2) < 0)
//						swap(all_boundary_V[0], all_boundary_V[1]);
//				}
//			}
//
//			//new_cell.site_normal /= index_of_adjacent_F.size();
//			//new_cell.site_normal.Normalized();
//			//new_cell.index_of_V = i;
//			new_cell.is_available = true;
//			new_cell.site = i;
//			new_cell.adjacent_cells = index_of_adjacent_V;	//�洢�ڽӶ��������
//			vector<Eigen::Vector3d> temp_lines;
//			all_lines.push_back(temp_lines);
//			for (int j = 0; j < all_boundary_V.size(); j++) {
//				Eigen::Vector3d temp_vec(all_boundary_V[j](0, 0), all_boundary_V[j](1, 0), all_boundary_V[j](2, 0));
//				all_lines[all_lines.size() - 1].push_back(temp_vec);
//				new_cell.all_points_in_polygon.push_back(Vector3(all_boundary_V[j](0, 0), all_boundary_V[j](1, 0), all_boundary_V[j](2, 0)));	//all_points_in_polygon�洢voronoi��߽綥��
//			}
//			all_voronoi_cells.push_back(new_cell);
//		}
//		else {
//			V_bottom.push_back(Vector3(V(i, 0), V(i, 1), V(i, 2)));
//			VoronoiCell new_cell;
//			new_cell.is_available = false;	//���ڸ߶���ֵ��Voronoi��Ԫ��Ϊ������
//			all_voronoi_cells.push_back(new_cell);
//		}
//	}
//	int end_time = clock();
//	std::cout << "&&&time&&& Sampling: " << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
//
//	Visual vis;
//	if (open_vis_voronoi == true)
//		vis.generateModelForRendering_8(all_lines, file_name + "_voronoi.obj");
//}
//
//
//void HybridManufacturing::GetVoronoiCells1()
//{
//	int start_time = clock();
//	double threshold_z = 2;
//	vector<vector<Eigen::Vector3d>> all_lines;
//	double min_z = MAX_D;
//	for (int i = 0; i < V.rows(); i++)
//		min_z = min(min_z, V(i, 2));
//	for (int i = 0; i < V.rows(); i++) {
//		if (V(i, 2) - min_z > threshold_z) {
//			VoronoiCell new_cell;
//			//new_cell.site_normal.m_x = new_cell.site_normal.m_y = new_cell.site_normal.m_z = 0;
//			vector<Eigen::MatrixXd> all_boundary_V;
//			vector<int> index_of_adjacent_F;
//			index_of_adjacent_F.clear();
//			vector<int> index_of_adjacent_V;
//			index_of_adjacent_V.clear();
//			for (int j = 0; j < F.rows(); j++) {
//				for (int k = 0; k < 3; k++) {
//					if (F(j, k) == i) {
//						index_of_adjacent_F.push_back(j);
//					}
//				}
//			}
//			int index_begin, index_end;
//			int index_first_begin;
//			for (int j = 0; j < 3; j++)
//				if (F(index_of_adjacent_F[0], j) == i) {
//					index_begin = F(index_of_adjacent_F[0], (j + 1) % 3);
//					index_end = F(index_of_adjacent_F[0], (j + 2) % 3);
//					index_first_begin = index_begin;
//					break;
//				}
//			index_of_adjacent_V.push_back(index_begin);
//			int cont_num = 1;
//			while (index_end != index_first_begin) {
//				for (int j = cont_num; j < index_of_adjacent_F.size(); j++) {
//					for (int k = 0; k < 3; k++) {
//						if (F(index_of_adjacent_F[j], k) == i) {
//							if (F(index_of_adjacent_F[j], (k + 1) % 3) == index_end) {
//								index_begin = F(index_of_adjacent_F[j], (k + 1) % 3);
//								index_end = F(index_of_adjacent_F[j], (k + 2) % 3);
//								index_of_adjacent_V.push_back(index_begin);
//								swap(index_of_adjacent_F[j], index_of_adjacent_F[cont_num]);
//								cont_num++;
//							}
//							else if (F(index_of_adjacent_F[j], (k + 2) % 3) == index_end) {
//								index_begin = F(index_of_adjacent_F[j], (k + 2) % 3);
//								index_end = F(index_of_adjacent_F[j], (k + 1) % 3);
//								index_of_adjacent_V.push_back(index_begin);
//								swap(index_of_adjacent_F[j], index_of_adjacent_F[cont_num]);
//								cont_num++;
//							}
//						}
//					}
//
//				}
//			}
//			for (int j = 0; j < index_of_adjacent_F.size(); j++) {
//				double na = (V(F(index_of_adjacent_F[j], 1), 1) - V(F(index_of_adjacent_F[j], 0), 1)) * (V(F(index_of_adjacent_F[j], 2), 2) - V(F(index_of_adjacent_F[j], 0), 2)) - (V(F(index_of_adjacent_F[j], 1), 2) - V(F(index_of_adjacent_F[j], 0), 2)) * (V(F(index_of_adjacent_F[j], 2), 1) - V(F(index_of_adjacent_F[j], 0), 1));
//				double nb = (V(F(index_of_adjacent_F[j], 1), 2) - V(F(index_of_adjacent_F[j], 0), 2)) * (V(F(index_of_adjacent_F[j], 2), 0) - V(F(index_of_adjacent_F[j], 0), 0)) - (V(F(index_of_adjacent_F[j], 1), 0) - V(F(index_of_adjacent_F[j], 0), 0)) * (V(F(index_of_adjacent_F[j], 2), 2) - V(F(index_of_adjacent_F[j], 0), 2));
//				double nc = (V(F(index_of_adjacent_F[j], 1), 0) - V(F(index_of_adjacent_F[j], 0), 0)) * (V(F(index_of_adjacent_F[j], 2), 1) - V(F(index_of_adjacent_F[j], 0), 1)) - (V(F(index_of_adjacent_F[j], 1), 1) - V(F(index_of_adjacent_F[j], 0), 1)) * (V(F(index_of_adjacent_F[j], 2), 0) - V(F(index_of_adjacent_F[j], 0), 0));
//				Eigen::Vector3d vn(na, nb, nc);
//				//new_cell.site_normal.m_x += vn[0];
//				//new_cell.site_normal.m_y += vn[1];
//				//new_cell.site_normal.m_z += vn[2];
//				Eigen::Matrix3d rotMatrix;
//
//				Eigen::Vector3d vectorBefore(0, 0, 1);
//				rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vn).toRotationMatrix();
//				vector<Eigen::MatrixXd> current_V;
//				Eigen::MatrixXd temp_V;
//				temp_V.resize(3, 1);
//				temp_V(0, 0) = V(i, 0);
//				temp_V(1, 0) = V(i, 1);
//				temp_V(2, 0) = V(i, 2);
//				current_V.push_back(temp_V);
//				temp_V(0, 0) = V(index_of_adjacent_V[j], 0);
//				temp_V(1, 0) = V(index_of_adjacent_V[j], 1);
//				temp_V(2, 0) = V(index_of_adjacent_V[j], 2);
//				current_V.push_back(temp_V);
//				temp_V(0, 0) = V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 0);
//				temp_V(1, 0) = V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 1);
//				temp_V(2, 0) = V(index_of_adjacent_V[(j + 1) % index_of_adjacent_V.size()], 2);
//				current_V.push_back(temp_V);
//				for (int k = 0; k < current_V.size(); k++) {
//					current_V[k] = rotMatrix.inverse() * current_V[k];
//				}
//
//				/*double x1 = current_V[0](0, 0), x2 = current_V[1](0, 0);
//				double y1 = current_V[0](1, 0), y2 = current_V[1](1, 0);
//				double B1 = 1, A1 = (x1 - x2) / (y1 - y2); double C1 = -((y1 + y2) / 2 + A1 * (x1 + x2) / 2);
//				x1 = current_V[0](0, 0), x2 = current_V[2](0, 0);
//				y1 = current_V[0](1, 0), y2 = current_V[2](1, 0);
//				double B2 = 1, A2 = (x1 - x2) / (y1 - y2); double C2 = -((y1 + y2) / 2 + A2 * (x1 + x2) / 2);
//
//				Eigen::MatrixXd new_V;
//				new_V.resize(3, 1);
//				new_V(0, 0) = (C2 * B1 - C1 * B2) / (A1 * B2 - A2 * B1);
//				new_V(1, 0) = (C1 * A2 - C2 * A1) / (A1 * B2 - A2 * B1);
//				new_V(2, 0) = current_V[0](2, 0);*/
//
//				Eigen::MatrixXd new_V;
//				new_V.resize(3, 1);
//				new_V(0, 0) = (current_V[0](0, 0) + current_V[1](0, 0) + current_V[2](0, 0)) / 3;
//				new_V(1, 0) = (current_V[0](1, 0) + current_V[1](1, 0) + current_V[2](1, 0)) / 3;
//				new_V(2, 0) = current_V[0](2, 0);
//
//				rotMatrix = Eigen::Quaterniond::FromTwoVectors(vn, vectorBefore).toRotationMatrix();
//				new_V = rotMatrix.inverse() * new_V;
//				all_boundary_V.push_back(new_V);
//
//				if (all_boundary_V.size() == 2) {
//					Eigen::Vector3d v1 = Eigen::Vector3d(V(i, 0), V(i, 1), V(i, 2));
//					Eigen::Vector3d v2 = Eigen::Vector3d(all_boundary_V[0](0, 0), all_boundary_V[0](1, 0), all_boundary_V[0](2, 0));
//					Eigen::Vector3d v3 = Eigen::Vector3d(all_boundary_V[1](0, 0), all_boundary_V[1](1, 0), all_boundary_V[1](2, 0));
//					//double ans = (v2.x() - v1.x()) * (v2.y() - v3.y()) - (v2.y() - v1.y()) * (v2.x() - v3.x());
//					//if (ans > 0)	//is clockwise
//					//	swap(v2, v3);
//					double na = (v2.y() - v1.y()) * (v3.z() - v1.z()) - (v2.z() - v1.z()) * (v3.y() - v1.y());
//					double nb = (v2.z() - v1.z()) * (v3.x() - v1.x()) - (v2.x() - v1.x()) * (v3.z() - v1.z());
//					double nc = (v2.x() - v1.x()) * (v3.y() - v1.y()) - (v2.y() - v1.y()) * (v3.x() - v1.x());
//					Eigen::Vector3d vn_2(na, nb, nc);
//					if (vn.dot(vn_2) < 0)
//						swap(all_boundary_V[0], all_boundary_V[1]);
//				}
//			}
//
//			//new_cell.site_normal /= index_of_adjacent_F.size();
//			//new_cell.site_normal.Normalized();
//			//new_cell.index_of_V = i;
//			new_cell.is_available = true;
//			new_cell.site = i;
//			new_cell.adjacent_cells = index_of_adjacent_V;
//			vector<Eigen::Vector3d> temp_lines;
//			all_lines.push_back(temp_lines);
//			for (int j = 0; j < all_boundary_V.size(); j++) {
//				Eigen::Vector3d temp_vec(all_boundary_V[j](0, 0), all_boundary_V[j](1, 0), all_boundary_V[j](2, 0));
//				all_lines[all_lines.size() - 1].push_back(temp_vec);
//				new_cell.all_points_in_polygon.push_back(Vector3(all_boundary_V[j](0, 0), all_boundary_V[j](1, 0), all_boundary_V[j](2, 0)));
//			}
//			all_voronoi_cells.push_back(new_cell);
//		}
//		else {
//			V_bottom.push_back(Vector3(V(i, 0), V(i, 1), V(i, 2)));
//			VoronoiCell new_cell;
//			new_cell.is_available = false;
//			all_voronoi_cells.push_back(new_cell);
//		}
//	}
//	int end_time = clock();
//	std::cout << "&&&time&&& Sampling: " << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
//
//	Visual vis;
//	if (open_vis_voronoi == true)
//		vis.generateModelForRendering_8(all_lines, file_name + "_voronoi.obj");
//}