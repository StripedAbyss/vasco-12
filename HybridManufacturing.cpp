#include "HybridManufacturing.h"
#include<stdlib.h>
#include<igl/readOBJ.h>
#include<igl/writeOBJ.h>
#include <array>

namespace
{
	inline void PrepareToolForCollision(cutter& tool)
	{
		tool.cylinder_height_threshold = tool.cylinder_height + tool.ball_r;
		tool.carriage_check_radius_sq = (tool.carriage_r + 5.0) * (tool.carriage_r + 5.0);
		tool.cylinder_check_radius_sq = (tool.cylinder_r + 5.0) * (tool.cylinder_r + 5.0);
		tool.cylinder_r_sq = tool.cylinder_r * tool.cylinder_r;
		tool.carriage_r_sq = tool.carriage_r * tool.carriage_r;
		tool.total_height = tool.cylinder_height + tool.ball_r + tool.carriage_height;
	}

	inline Eigen::Vector3d ToVector3(const Eigen::MatrixXd& vec)
	{
		return { vec(0, 0), vec(1, 0), vec(2, 0) };
	}

	inline Eigen::Vector3d ComputeToolCenter(const Eigen::Vector3d& point, const Eigen::Vector3d& normal, double radius)
	{
		return point + radius * normal;
	}
}

HybridManufacturing::HybridManufacturing(std::string file_name, std::string suf, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd N)
{
	this->V = V;
	this->F = F;
	this->Normals = N;
	this->file_name = file_name;
	this->suf = suf;
	InitializePolyscope();
}

HybridManufacturing::~HybridManufacturing()
{
}

std::vector<TRiangle> HybridManufacturing::FilterSurfaceRemoveTriangles(
	const Slicer_2& slicer,
	const std::vector<TRiangle>& remove_triangles) const
{
	std::vector<TRiangle> surface_triangles;
	surface_triangles.reserve(remove_triangles.size());

	const size_t face_cnt = Normals.rows();
	for (size_t i = 0; i < remove_triangles.size(); i++) {
		Eigen::Vector3d v1, v2, v3;
		v1.x() = slicer.positions[remove_triangles[i][0]][0];
		v1.y() = slicer.positions[remove_triangles[i][0]][1];
		v1.z() = slicer.positions[remove_triangles[i][0]][2];
		v2.x() = slicer.positions[remove_triangles[i][1]][0];
		v2.y() = slicer.positions[remove_triangles[i][1]][1];
		v2.z() = slicer.positions[remove_triangles[i][1]][2];
		v3.x() = slicer.positions[remove_triangles[i][2]][0];
		v3.y() = slicer.positions[remove_triangles[i][2]][1];
		v3.z() = slicer.positions[remove_triangles[i][2]][2];
		double ans = (v2.x() - v1.x()) * (v2.y() - v3.y()) - (v2.y() - v1.y()) * (v2.x() - v3.x());
		if (ans > 0) {
			swap(v2, v3);
		}
		double na = (v2.y() - v1.y()) * (v3.z() - v1.z()) - (v2.z() - v1.z()) * (v3.y() - v1.y());
		double nb = (v2.z() - v1.z()) * (v3.x() - v1.x()) - (v2.x() - v1.x()) * (v3.z() - v1.z());
		double nc = (v2.x() - v1.x()) * (v3.y() - v1.y()) - (v2.y() - v1.y()) * (v3.x() - v1.x());

		Eigen::Vector3d normal_vector(na, nb, nc);
		normal_vector.normalize();

		bool jud_surface = false;
		for (size_t j = 0; j < face_cnt; j++) {
			double n_dot = normal_vector.dot(Normals.row(j));
			if (fabs(n_dot) < 0.999) {
				continue;
			}
			Eigen::Vector3d v1_s(v1.x(), v1.y(), v1.z());
			Eigen::Vector3d v2_s(v2.x(), v2.y(), v2.z());
			Eigen::Vector3d v3_s(v3.x(), v3.y(), v3.z());

			bool b1 = vasco::core::isPointInTriangle(v1_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));
			bool b2 = vasco::core::isPointInTriangle(v2_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));
			bool b3 = vasco::core::isPointInTriangle(v3_s, V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)));

			if (b1 && b2 && b3) {
				jud_surface = true;
				break;
			}
		}

		if (jud_surface) {
			surface_triangles.push_back(remove_triangles[i]);
		}
	}

	return surface_triangles;
}

void HybridManufacturing::InitializeVoronoi(const std::vector<VoronoiCell>& cells,
	const std::vector<Eigen::Vector3d>& bottom)
{
	all_voronoi_cells = cells;
	V_bottom = bottom;
}

void HybridManufacturing::InitializePolyscope()
{
	polyscope::init();
	polyscope::view::setUpDir(polyscope::UpDir::ZUp);
}

int HybridManufacturing::CollisionDetectionForSubtractiveManufacturing(cutter the_nozzle)
{
	clock_t start_time, end_time;
	clock_t start_time_2, end_time_2;
	clock_t start_time_3, end_time_3;
	clock_t start_time_4, end_time_4;
	double sum_time = 0;
	int test_cont_num = 0;

	start_time = clock();
	std::cout << "Doing collision detection for subtractive manufacturing......" << endl;
	string file_name_2 = file_name;
	int cont_num = 0;
	sampling_subtractive.OrientationSamplePoints();	//sampling_subtractive生成球面采样点
	Eigen::Matrix3d rotMatrix;
	std::vector<bool> flag_accessible_points;	//true值代表对应的voronoi单元在某个采样方向下可达
	std::vector<bool> flag_covering_points;	//true值代表对应的voronoi单元在某个采样方向下被喷嘴覆盖
	//double nozzle_par = (the_nozzle.upper_surface_r - the_nozzle.lowwer_surface_r) / the_nozzle.nozzle_H_half;

	flag_accessible_points.resize(all_voronoi_cells.size());
	flag_covering_points.resize(all_voronoi_cells.size());
	flag_voronoi_has_been_printed.resize(all_voronoi_cells.size());
	for (int i = 0; i < all_voronoi_cells.size(); i++) {
		flag_accessible_points[i] = false;
		flag_voronoi_has_been_printed[i] = true;
		flag_covering_points[i] = false;
	}

	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {	//枚举所有采样方向

		///////////////////rotate/////////////////////
		std::vector<Eigen::MatrixXd> temp_V;
		temp_V.resize(V.rows());	//temp_V存储旋转后的模型顶点
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i].resize(3, 1);
			temp_V[i](0, 0) = V.row(i).x();
			temp_V[i](1, 0) = V.row(i).y();
			temp_V[i](2, 0) = V.row(i).z();
		}
		vector<std::vector<Eigen::MatrixXd>> temp_new_V;	//temp_new_V存储旋转后的voronoi面边界顶点
		temp_new_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_new_V[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
			for (int k = 0; k < temp_new_V[i].size(); k++) {
				temp_new_V[i][k].resize(3, 1);
				temp_new_V[i][k](0, 0) = all_voronoi_cells[i].all_points_in_polygon[k].x();
				temp_new_V[i][k](1, 0) = all_voronoi_cells[i].all_points_in_polygon[k].y();
				temp_new_V[i][k](2, 0) = all_voronoi_cells[i].all_points_in_polygon[k].z();
			}
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();	//计算从z轴到采样方向的旋转矩阵
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix.inverse() * temp_V[i];
		for (int i = 0; i < V.rows(); i++)
			for (int j = 0; j < temp_new_V[i].size(); j++)
				temp_new_V[i][j] = rotMatrix.inverse() * temp_new_V[i][j];
		//////////////////////////////////////////////


		/////////////////calculate normal//////////////////
		//vector<vector<Vector3>> all_normal_of_triangles_in_cells;
		vector<Eigen::Vector3d> all_normal_of_cells;	//存储每个voronoi单元的法向量
		//all_normal_of_triangles_in_cells.resize(all_voronoi_cells.size());
		all_normal_of_cells.resize(all_voronoi_cells.size());

		vector<vector<Eigen::Vector3d>> temp_vis;
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			if (all_voronoi_cells[i].is_available == true) {
				//all_normal_of_triangles_in_cells[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
				all_normal_of_cells[i].x() = all_normal_of_cells[i].y() = all_normal_of_cells[i].z() = 0;
				for (int j = 0; j < 1; j++) {	//计算voronoi单元法向量时只用第一个三角形？v1是site，v2、v3是边界顶点
					Eigen::Vector3d v1 = ToVector3(temp_V[all_voronoi_cells[i].site]);
					Eigen::Vector3d v2 = ToVector3(temp_new_V[i][j]);
					Eigen::Vector3d v3 = ToVector3(temp_new_V[i][(j + 1) % all_voronoi_cells[i].all_points_in_polygon.size()]);
					Eigen::Vector3d vn = (v2 - v1).cross(v3 - v1);
					//all_normal_of_triangles_in_cells[i][j] = vn;
					all_normal_of_cells[i].x() += vn.x();
					all_normal_of_cells[i].y() += vn.y();
					all_normal_of_cells[i].z() += vn.z();
				}
				all_normal_of_cells[i] /= all_voronoi_cells[i].all_points_in_polygon.size();
				all_normal_of_cells[i].normalize();	//voronoi单元法向量归一化

				all_voronoi_cells[i].all_normal_in_all_ori.push_back(all_normal_of_cells[i]);	//将该采样方向下的法向量存入all_voronoi_cells的all_normal_in_all_ori中

				vector<Eigen::Vector3d> temp_vec;
				temp_vis.push_back(temp_vec);
				Eigen::Vector3d v_site(
					temp_V[all_voronoi_cells[i].site](0, 0),
					temp_V[all_voronoi_cells[i].site](1, 0),
					temp_V[all_voronoi_cells[i].site](2, 0));
				Eigen::Vector3d v_normal = v_site + all_normal_of_cells[i] * 3.0;
				temp_vis[temp_vis.size() - 1].push_back(v_site);
				temp_vis[temp_vis.size() - 1].push_back(v_normal);
			}
		}
		if (ori == 0) {
			Visual vis_normal;//存储法向量的可视化
			//vis_normal.generateModelForRendering_9(temp_vis, file_name);

		}
		//////////////////////////////////////////////////


		///////////////////collision detection////////////////////////
		int cont_accessible_points = 0, cont_unaccessible_points = 0;

		vector<double> max_z_of_cells(all_voronoi_cells.size());	//max_z_of_cells存储每个voronoi单元边界顶点(all_points_in_polygon)的最大z值
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			max_z_of_cells[i] = MIN_D;
			for (int j = 0; j < all_voronoi_cells[i].all_points_in_polygon.size(); j++)
				max_z_of_cells[i] = max(max_z_of_cells[i], temp_new_V[i][j](2, 0));
		}

		PrepareToolForCollision(the_nozzle);

		for (int i = 0; i < all_voronoi_cells.size(); i++) {	//枚举所有voronoi单元，判断该单元在当前采样方向ori下是否可达
			if (flag_accessible_points[i] == true) {	//如果该单元在之前某个采样方向下已经可达，则跳过
				cont_accessible_points++;
				continue;
			}
			cont_unaccessible_points++;

			if (!all_voronoi_cells[i].is_available) {
				continue;
			}

			// 计算刀尖球中心点坐标
			Eigen::Vector3d center_point = ComputeToolCenter(
				ToVector3(temp_V[all_voronoi_cells[i].site]),
				all_normal_of_cells[i],
				the_nozzle.cylinder_r);

			bool jud_collision = false;
			for (int ii = 0; ii < all_voronoi_cells.size(); ii++) {
				if (i == ii || !all_voronoi_cells[ii].is_available) {
					continue;
				}

				if (CheckToolCollisionWithCell(center_point, temp_new_V[ii], max_z_of_cells[ii], the_nozzle)) {
					jud_collision = true;
					break;
				}
			}

			if (!jud_collision) {
				flag_accessible_points[i] = true;
			}
		}
		num_inaccessible_points = cont_unaccessible_points;
		std::cout << "id of orientation:" << ori << endl;
		std::cout << "number of accessible points:" << cont_accessible_points << endl;
		std::cout << "number of unaccessible points:" << cont_unaccessible_points << endl << endl;
		//break;
	}
	end_time = clock();


	//*******************************************//
	/////////////////find area S//////////////////
	//*******************************************//
	start_time_2 = clock();

	int cont_number_2 = 0;
	std::vector<Eigen::MatrixXd> vis_red_points;	//vis_red_points存储不可达voronoi单元的边界顶点坐标，用于可视化
	std::vector<Eigen::MatrixXd> vis_green_points;	//vis_green_points存储可达voronoi单元的边界顶点坐标，用于可视化
	std::vector<Eigen::MatrixXd> temp_V_vis;	//temp_V_vis存储未旋转的模型顶点坐标，用于可视化
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {	//ori枚举所有采样方向
		std::vector<Eigen::MatrixXd> temp_V;	//temp_V存储旋转后的模型顶点坐标
		temp_V.resize(V.rows());
		temp_V_vis.resize(V.rows());	//temp_V_vis存储未旋转的模型顶点坐标
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i].resize(3, 1);
			temp_V[i](0, 0) = V.row(i).x();
			temp_V[i](1, 0) = V.row(i).y();
			temp_V[i](2, 0) = V.row(i).z();
		}
		vector<std::vector<Eigen::MatrixXd>> temp_new_V;	//temp_new_V存储旋转后的voronoi面边界顶点坐标
		temp_new_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_new_V[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
			for (int k = 0; k < temp_new_V[i].size(); k++) {
				temp_new_V[i][k].resize(3, 1);
				temp_new_V[i][k](0, 0) = all_voronoi_cells[i].all_points_in_polygon[k].x();
				temp_new_V[i][k](1, 0) = all_voronoi_cells[i].all_points_in_polygon[k].y();
				temp_new_V[i][k](2, 0) = all_voronoi_cells[i].all_points_in_polygon[k].z();
			}
		}
		temp_V_vis = temp_V;

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix.inverse() * temp_V[i];	//将模型顶点旋转到以采样方向为z轴的坐标系下
		for (int i = 0; i < V.rows(); i++)
			for (int j = 0; j < temp_new_V[i].size(); j++)
				temp_new_V[i][j] = rotMatrix.inverse() * temp_new_V[i][j];	//将voronoi面边界顶点旋转到以采样方向为z轴的坐标系下

		//cout << "b" << endl;
		/////////////////calculate normal//////////////////
		//vector<vector<Vector3>> all_normal_of_triangles_in_cells;
		vector<Eigen::Vector3d> all_normal_of_cells;	//仅计算不可达voronoi单元的法向量
		//all_normal_of_triangles_in_cells.resize(all_voronoi_cells.size());
		all_normal_of_cells.resize(all_voronoi_cells.size());
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			if (all_voronoi_cells[i].is_available == true && flag_accessible_points[i] == false) {	//仅计算不可达voronoi单元的法向量
				//all_normal_of_triangles_in_cells[i].resize(all_voronoi_cells[i].all_points_in_polygon.size());
				all_normal_of_cells[i].x() = all_normal_of_cells[i].y() = all_normal_of_cells[i].z() = 0;
				for (int j = 0; j < 1; j++) {
					Eigen::Vector3d v1 = ToVector3(temp_V[all_voronoi_cells[i].site]);
					Eigen::Vector3d v2 = ToVector3(temp_new_V[i][j]);
					Eigen::Vector3d v3 = ToVector3(temp_new_V[i][(j + 1) % all_voronoi_cells[i].all_points_in_polygon.size()]);
					Eigen::Vector3d vn = (v2 - v1).cross(v3 - v1);
					//all_normal_of_triangles_in_cells[i][j] = vn;
					all_normal_of_cells[i].x() += vn.x();
					all_normal_of_cells[i].y() += vn.y();
					all_normal_of_cells[i].z() += vn.z();
				}
				all_normal_of_cells[i] /= all_voronoi_cells[i].all_points_in_polygon.size();
				all_normal_of_cells[i].normalize();
			}
		}
		//////////////////////////////////////////////////
		//cout << "c" << endl;
		//ofstream vis_unaccessible_points(file_name + "_unaccessible_points.obj");

		file_name_2 = file_name + to_string(cont_num);
		vector<double> max_z_of_cells(all_voronoi_cells.size());	//max_z_of_cells存储每个voronoi单元边界顶点(all_points_in_polygon)的最大z值
		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			max_z_of_cells[i] = MIN_D;
			for (int j = 0; j < all_voronoi_cells[i].all_points_in_polygon.size(); j++)
				max_z_of_cells[i] = max(max_z_of_cells[i], temp_new_V[i][j](2, 0));
		}

		int cont_unaccessible_points = 0;

		PrepareToolForCollision(the_nozzle);

		for (int i = 0; i < all_voronoi_cells.size(); i++) {
			if (flag_accessible_points[i] == true) {	//跳过可达voronoi单元
				continue;
			}

			if (!all_voronoi_cells[i].is_available) {
				continue;
			}

			// 计算不可达voronoi单元i的刀尖球中心点坐标
			Eigen::Vector3d center_point = ComputeToolCenter(
				ToVector3(temp_V[all_voronoi_cells[i].site]),
				all_normal_of_cells[i],
				the_nozzle.cylinder_r);

			int index_insert = 0;
			if (ori == 0) {	//在首个方向 ori==0 时，初始化该不可达单元的 area_S 容器与索引映射
				vector<area_S> temp_vec_area_s;
				//temp_vec_area_s.resize(1);
				all_the_area_S.push_back(temp_vec_area_s);
				map_S_and_vertex.insert({ cont_unaccessible_points,i });
				map_S_and_vertex_inv.insert({ i, cont_unaccessible_points });
				cont_unaccessible_points++;	//记录不可达单元数量,第(cont_unaccessible_points-1)个不可达voronoi单元的索引为i
			}
			else {
				index_insert = map_S_and_vertex_inv[i];
			}

			//vis_unaccessible_points << "v " << temp_V[i](0,0) << " " << temp_V[i](1, 0) << " " << temp_V[i](2, 0) << endl;
			start_time_4 = clock();
			if (ori == 0)
				vis_red_points.push_back(temp_V_vis[i]);	//将不可达voronoi单元i的网格顶点存入vis_red_points，用于可视化
			end_time_4 = clock();
			sum_time += (double(end_time_4) - double(start_time_4)) / CLOCKS_PER_SEC;

			for (int ii = 0; ii < all_voronoi_cells.size(); ii++) {	//枚举所有voronoi单元，判断刀尖放在当前单元i、刀具方向是ori方向时，是否与其他单元发生碰撞
				if (i == ii || !all_voronoi_cells[ii].is_available) {
					continue;
				}

				if (CheckToolCollisionWithCell(center_point, temp_new_V[ii], max_z_of_cells[ii], the_nozzle)) {
					area_S temp_area_S(ii, ori);
					if (ori == 0) {
						all_the_area_S.back().push_back(temp_area_S);
					}
					else {
						all_the_area_S[index_insert].push_back(temp_area_S);
						cont_number_2++;
					}
					test_cont_num++;
				}
			}
		}

		//cout << "d" << endl;
		Visual Vis;
		//Vis.generateModelForRendering(temp_layers, file_name_2);
		//Vis.generateModelForRendering_3(vectorAfter, file_name_2, vis_points);
		//break;

		cont_num++;
	}
	//*******************************************//
	//cout << endl << all_the_area_S.size() << endl;
	//cout << endl << cont_number_2 << endl;
	end_time_2 = clock();

	////////////get covering points(green points)////////////
	//only consider covering for area S,save the index of all_the_area_S
	start_time_3 = clock();
	cout << test_cont_num << "&&&&" << endl;
	int cont_covering_points = 0;
	for (int i = 0; i < all_the_area_S.size(); i++) {	//枚举第i个不可达voronoi单元的area_S
		for (int j = 0; j < all_the_area_S[i].size(); j++) {	//枚举里面的所有area_S元素
			if (flag_covering_points[all_the_area_S[i][j].pointId] == false) {	//如果该voronoi单元还没有在flag_covering_points中被标记
				map_covering_points_and_vertex.insert({ cont_covering_points,all_the_area_S[i][j].pointId });	//建立被碰撞点与voronoi单元索引的映射
				map_covering_points_and_vertex_inv.insert({ all_the_area_S[i][j].pointId, cont_covering_points });
				vector<area_S> temp_vec_area_s;
				all_the_covering_points.push_back(temp_vec_area_s);	//在all_the_covering_points中添加一个新的容器

				area_S temp_covering_point(i, all_the_area_S[i][j].oriId);	//
				all_the_covering_points[all_the_covering_points.size() - 1].push_back(temp_covering_point);
				flag_covering_points[all_the_area_S[i][j].pointId] = true;
				cont_covering_points++;
			}
			else {
				int index_insert = map_covering_points_and_vertex_inv[all_the_area_S[i][j].pointId];
				area_S temp_covering_point(i, all_the_area_S[i][j].oriId);
				all_the_covering_points[index_insert].push_back(temp_covering_point);
			}
		}
	}
	int max_size = -MAX_I;
	for (int i = 0; i < all_the_covering_points.size(); i++) {
		max_size = max(max_size, int(all_the_covering_points[i].size()));
		vis_green_points.push_back(temp_V_vis[map_covering_points_and_vertex[i]]);	//将覆盖点对应的voronoi单元顶点存入vis_green_points，用于可视化
	}
	vector<double> color_map;
	color_map.resize(all_the_covering_points.size());
	for (int i = 0; i < all_the_covering_points.size(); i++)
		color_map[i] = double(all_the_covering_points[i].size()) / double(max_size);	//归一化颜色值，便于可视化，被更多个area S覆盖的点颜色更亮
	end_time_3 = clock();
	/*vector<vector<area_S>> temp_all_the_covering_points = all_the_covering_points;
	map<int, int> temp_map_covering_points_and_vertex = map_covering_points_and_vertex;
	for(int i = 0;i< temp_all_the_covering_points.size();i++)
		for (int j = i + 1; j < temp_all_the_covering_points.size(); j++)
			if (temp_all_the_covering_points[i].size() < temp_all_the_covering_points[j].size()) {
				swap(temp_all_the_covering_points[i], temp_all_the_covering_points[j]);
				swap(temp_map_covering_points_and_vertex[i], temp_map_covering_points_and_vertex[j]);
			}

	for(int i =0;i < temp_all_the_covering_points.size() / 10;i++)
		vis_green_points.push_back(temp_V_vis[temp_map_covering_points_and_vertex[i]]);*/
		/////////////////////////////////////////////////////////
	if (open_vis_red_points == true)
		createRedBalls(file_name, vis_red_points);
	if (open_vis_green_points == true)
		createGreenBalls(file_name, vis_green_points, color_map);

	std::cout << "Collision Detection For Subtractive Manufacturing have done" << endl;
	/*std::cout << "&&&time&&& Collision detection: " << double(end_time - start_time) / CLOCKS_PER_SEC << std::endl;
	std::cout << "&&&time&&& Subtractive dependency graph: " << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << std::endl;
	std::cout << "&&&time&&& Green points generation: " << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << std::endl;*/
	time_build_subtractive_graph = double(end_time - start_time) / CLOCKS_PER_SEC + double(end_time_2 - start_time_2) / CLOCKS_PER_SEC + double(end_time_3 - start_time_3) / CLOCKS_PER_SEC;
	std::cout << "&&&time&&& Build subtractive graph: " << double(end_time - start_time) / CLOCKS_PER_SEC + double(end_time_2 - start_time_2) / CLOCKS_PER_SEC + double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << std::endl;

	if (all_the_area_S.size() == 0) {
		cout << file_name << " 无不可达点！！！" << endl;

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXi N;

		igl::readOBJ(file_name + ".obj", V, F);
		igl::writeOBJ(file_name + "-1_0_current.obj", V, F);
		return -1;
	}
	return 0;
}

void HybridManufacturing::GetALLFragileVertex(SAMPLE_ON_BALL sampling)
{
	for (int ori = 0; ori < sampling.sample_points.size(); ori++) {
		//rotating the blocks and then slicing//
		vector<Eigen::MatrixXd> temp_V;
		temp_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			temp_V[i].resize(3, 1);
			temp_V[i](0, 0) = V.row(i).x();
			temp_V[i](1, 0) = V.row(i).y();
			temp_V[i](2, 0) = V.row(i).z();
		}
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling.sample_points[ori]);
		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < V.rows(); i++)
			temp_V[i] = rotMatrix.inverse() * temp_V[i];

		vector<Eigen::Vector3d> normal_V;
		vector<bool> temp_fragile_V;
		normal_V.resize(V.rows());
		temp_fragile_V.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			normal_V[i].x() = normal_V[i].y() = normal_V[i].z() = 0;
		}
		for (int i = 0; i < F.rows(); i++) {
			VEctor v1, v2, v3;
			v1[0] = temp_V[F(i, 0)](0, 0);
			v1[1] = temp_V[F(i, 0)](1, 0);
			v1[2] = temp_V[F(i, 0)](2, 0);
			v2[0] = temp_V[F(i, 1)](0, 0);
			v2[1] = temp_V[F(i, 1)](1, 0);
			v2[2] = temp_V[F(i, 1)](2, 0);
			v3[0] = temp_V[F(i, 2)](0, 0);
			v3[1] = temp_V[F(i, 2)](1, 0);
			v3[2] = temp_V[F(i, 2)](2, 0);
			Eigen::Vector3d temp_normal;
			temp_normal = faceNormal(v1, v2, v3);
			temp_normal.normalize();
			normal_V[F(i, 0)] += temp_normal;
			normal_V[F(i, 1)] += temp_normal;
			normal_V[F(i, 2)] += temp_normal;
		}
		for (int i = 0; i < V.rows(); i++) {
			normal_V[i].normalize();
			if (normal_V[i].z() > 0.95)
				temp_fragile_V[i] = true;
			else
				temp_fragile_V[i] = false;
		}
		is_fragile_V.push_back(temp_fragile_V);
	}
}

void HybridManufacturing::detect_collision_with_printing_platform(
	int& index,
	vector<int>& candidate_nodes,
	OrientationScores& all_calculated_value,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d ori_now,
	nozzle the_nozzle)
{

	Eigen::Vector3d center_point(0.0, 0.0, 0.0);
	for (int i = 0; i < V_bottom.size(); i++) {
		center_point += V_bottom[i];
	}
	center_point /= V_bottom.size();
	double circle_r = 160;

	if (ori_now.x() == 0 && ori_now.y() == 0 && ori_now.z() == 1)
		return;
	//if (all_calculated_value[index].number_of_remaining_face < terminate_threshold_of_number_of_faces) 
		//return;
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(ori_now, vectorBefore).toRotationMatrix();
	vector<vector<Eigen::MatrixXd>> rotate_all_cut_layers;
	rotate_all_cut_layers.resize(all_cut_layers.size());
	for (int i = 0; i < all_cut_layers.size(); i++) {
		rotate_all_cut_layers[i].resize(all_cut_layers[i].size());
		for (int j = 0; j < all_cut_layers[i].size(); j++) {
			rotate_all_cut_layers[i][j].resize(3, 1);
			rotate_all_cut_layers[i][j](0, 0) = all_cut_layers[i][j].x;
			rotate_all_cut_layers[i][j](1, 0) = all_cut_layers[i][j].y;
			rotate_all_cut_layers[i][j](2, 0) = all_cut_layers[i][j].z;
			rotate_all_cut_layers[i][j] = rotMatrix.inverse() * rotate_all_cut_layers[i][j];
		}
	}

	for (int i = 0; i < rotate_all_cut_layers.size(); i++) {
		for (int j = 0; j < rotate_all_cut_layers[i].size(); j++) {
			if ((abs(rotate_all_cut_layers[i][j](2, 0) - center_point.z()) < the_nozzle.lowwer_surface_r) && (pow(rotate_all_cut_layers[i][j](0, 0) - center_point.x(), 2) + pow(rotate_all_cut_layers[i][j](1, 0) - center_point.y(), 2) - pow(circle_r, 2) < 0)) {
				candidate_nodes.erase(candidate_nodes.begin() + index);
				all_calculated_value.erase(all_calculated_value.begin() + index);
				index--;
				return;
			}
		}
	}
}

all_value HybridManufacturing::GainMesh(
	Slicer_2& slicer,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d vector_after,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	vector<int> all_cut_layers_dependency_layer,
	bool flag_is_continue_block,
	int id_continue)
{
	clock_t start_time, end_time;

	double sum_area = 0;
	double sum_value_of_manufacturing_require = 0;
	double sum_not_self_support_area = 0;
	//double area_threshold = 5.0;
	double ratio_of_area_threshold = 0.05;  //0.015 //0.15  //0.08   //0.03  //0.05
	double max_self_support_slope_angle;
	//cout << "#####" << height_of_beam_search << " " << id_continue << endl;
	//if(height_of_beam_search != 7 || id_continue != 1)s
	//	max_self_support_slope_angle = PI / 10;  //PI/3.6  
	//else
	max_self_support_slope_angle = PI / 3.6;  //PI/3.6  
	double sum_area_threshold = 1000;  //200 //1000
	if (vector_after.x() == 0 && vector_after.y() == 0) {
		sum_area_threshold = 100000;
		ratio_of_area_threshold = 0.1;
	}

	/*if(height_of_beam_search == 4)
		ratio_of_area_threshold = 0.2;*/

	std::vector<bool> jud_triangle_have_been_added;
	vector<vector<int>> cutting_plane_points;
	cutting_plane_points.resize(all_cut_layers.size());
	//////sort cut layers//////
	for (int i = 0; i < all_cut_layers.size(); i++)
		for (int j = i + 1; j < all_cut_layers.size(); j++) {

			if (all_cut_layers[i][0].z > all_cut_layers[j][0].z) {
				swap(all_cut_layers[i], all_cut_layers[j]);
				swap(all_cut_layers_dependency_layer[i], all_cut_layers_dependency_layer[j]);
			}
		}
	///////////////////////////


	//need rotate first//
	Eigen::Vector3d vectorBefore(0, 0, 1);
	Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vector_after).toRotationMatrix();
	for (int i = 0; i < slicer.positions.size(); i++) {
		Eigen::MatrixXd temp_V;
		temp_V.resize(3, 1);
		temp_V(0, 0) = slicer.positions[i][0];
		temp_V(1, 0) = slicer.positions[i][1];
		temp_V(2, 0) = slicer.positions[i][2];
		temp_V = rotMatrix.inverse() * temp_V;
		slicer.positions[i][0] = temp_V(0, 0);
		slicer.positions[i][1] = temp_V(1, 0);
		slicer.positions[i][2] = temp_V(2, 0);
	}
	/////////////////////


	//--------------------cut-----------------------//
	slicer.normal[0] = 0;
	slicer.normal[1] = 0;
	slicer.normal[2] = 1;
	//cout << endl << "pp" << endl;
	for (int i = 0; i < all_cut_layers.size(); i++) {
		slicer.origin[0] = 0;
		slicer.origin[1] = 0;
		slicer.origin[2] = all_cut_layers[i][0].z;
		slicer.cut();
	}
	//cout << endl << "qq" << endl;
	std::vector<TRiangle> ori_triangle = slicer.triangles;
	jud_triangle_have_been_added.resize(slicer.triangles.size());
	for (int i = 0; i < slicer.triangles.size(); i++)
		jud_triangle_have_been_added[i] = false;

	Slicer_2 all_slicer;

	all_slicer.positions = slicer.positions;
	all_slicer.triangles = ori_triangle;


	//--------------------save candidate_triangles-----------------------//
	int current_index = 0;
	std::vector<VEctor> min_z_point;
	std::vector<double> min_z_triangle;
	vector<int> index_of_min_point_in_triangle;
	std::vector<TRiangle> candidate_triangles;
	std::vector<int> id_triangles;
	std::vector<TRiangle> remove_triangles;
	candidate_triangles.clear();
	id_triangles.clear();
	double min_z_all_cut_layers = 9999999;
	for (int t = 0; t < all_cut_layers.size(); t++)
		min_z_all_cut_layers = min(min_z_all_cut_layers, all_cut_layers[t][0].z);
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;
		for (int k = 0; k < 3; k++) {
			temp_min_z_triangle = min(all_slicer.positions[all_slicer.triangles[i][k]][2], temp_min_z_triangle);
		}
		if (temp_min_z_triangle + 0.001 >= min_z_all_cut_layers) {
			candidate_triangles.push_back(all_slicer.triangles[i]);
			id_triangles.push_back(i);
		}
	}
	for (int i = 0; i < candidate_triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;
		VEctor temp_min_z_point = all_slicer.positions[candidate_triangles[i][0]];
		int temp_index_of_min_point_in_triangle;
		for (int k = 0; k < 3; k++) {
			if (all_slicer.positions[candidate_triangles[i][k]][2] < temp_min_z_triangle) {
				temp_min_z_triangle = all_slicer.positions[candidate_triangles[i][k]][2];
				temp_min_z_point = all_slicer.positions[candidate_triangles[i][k]];
				temp_index_of_min_point_in_triangle = candidate_triangles[i][k];
			}
		}
		min_z_triangle.push_back(temp_min_z_triangle);
		min_z_point.push_back(temp_min_z_point);
		index_of_min_point_in_triangle.push_back(temp_index_of_min_point_in_triangle);
	}

	quick_sort(candidate_triangles, 0, candidate_triangles.size() - 1, id_triangles, min_z_triangle, min_z_point);


	////--------------------save OPP_triangles one by one-----------------------//
	vector<vector<TRiangle>> all_furcation_of_blocks;
	all_furcation_of_blocks.resize(all_cut_layers.size());
	vector<int> id_of_furcation_of_blocks;
	vector<int> save_current_index;
	for (int t = 0; t < all_cut_layers.size(); t++) {
		current_index = 0;   //该方式也许比较慢
		double boundary_bottom = 999999, boundary_left = 999999, boundary_top = -999999, boundary_right = -999999;
		for (int i = 0; i < all_cut_layers[t].size(); i++) {
			boundary_top = std::max(boundary_top, all_cut_layers[t][i].y);
			boundary_bottom = std::min(boundary_bottom, all_cut_layers[t][i].y);
			boundary_right = std::max(boundary_right, all_cut_layers[t][i].x);
			boundary_left = std::min(boundary_left, all_cut_layers[t][i].x);
		}
		Eigen::Vector2d current_triangle_point;
		Eigen::Vector2d current_layer_point;
		for (; current_index < candidate_triangles.size(); current_index++) {
			if (abs(min_z_triangle[current_index] - all_cut_layers[t][0].z) > 0.001) {
				if (min_z_triangle[current_index] > all_cut_layers[t][0].z)
					break;
			}
			else {
				bool jud_is_boundary_point = false;
				//for (int i = 0; i < all_cut_layers[t].size(); i++) {
				int cont_inside_boundary = 0;
				for (int k = 0; k < 3; k++) {
					if (all_slicer.positions[candidate_triangles[current_index][k]][0] + 0.01 >= boundary_left && all_slicer.positions[candidate_triangles[current_index][k]][0] - 0.01 <= boundary_right
						&& all_slicer.positions[candidate_triangles[current_index][k]][1] + 0.01 >= boundary_bottom && all_slicer.positions[candidate_triangles[current_index][k]][1] - 0.01 <= boundary_top) {
						cont_inside_boundary++;
					}
				}
				if (cont_inside_boundary >= 2)
					jud_is_boundary_point = true;
				/*if (min_z_point[current_index][0] + 0.01 >= boundary_left && min_z_point[current_index][0] - 0.01 <= boundary_right &&
					min_z_point[current_index][1] + 0.01 >= boundary_bottom && min_z_point[current_index][1] - 0.01 <= boundary_top) {
					jud_is_boundary_point = true;
					break;
				}*/
				//}
				if (jud_is_boundary_point == true) {
					for (int j = 0; j < all_cut_layers[t].size(); j++) {
						current_layer_point.x() = all_cut_layers[t][j].x;
						current_layer_point.y() = all_cut_layers[t][j].y;
						current_triangle_point.x() = min_z_point[current_index][0];
						current_triangle_point.y() = min_z_point[current_index][1];
						if ((current_layer_point - current_triangle_point).norm() < 0.2) {  //4.0
							//if (height_of_beam_search != 2)
							remove_triangles.push_back(candidate_triangles[current_index]);
							id_of_furcation_of_blocks.push_back(t);
							all_furcation_of_blocks[t].push_back(candidate_triangles[current_index]);

							Eigen::Vector3d face_normal = faceNormal(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							face_normal.normalize();
							Eigen::Vector3d base_normal(0, 0, 1);
							int jud_self_support = (face_normal.dot(base_normal) + sin(max_self_support_slope_angle) >= 0);
							double current_area = triangleArea(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							sum_area += current_area;
							//sum_value_of_manufacturing_require += current_area * calculate_manufacturing_require(all_slicer, remove_triangles);
							if (jud_self_support == false) {
								sum_not_self_support_area += current_area;
								if (sum_not_self_support_area / sum_area > ratio_of_area_threshold && sum_area >= sum_area_threshold) {
									all_value temp_all_value;
									temp_all_value.value_of_self_support = 0;
									return temp_all_value;
								}
							}

							save_current_index.push_back(current_index);
							jud_triangle_have_been_added[id_triangles[current_index]] = true;
							break;
						}
					}
				}
			}

		}
	}
	current_index = 0;  //current_index = 0;
	//std::cout << "**" << current_index << " " << candidate_triangles.size() << " " << remove_triangles.size() << endl;

	start_time = clock();
	while (1) {
		bool flag_break = false;
		for (int i = 0; i < remove_triangles.size(); i++)
		{
			if (jud_triangle_have_been_added[id_triangles[current_index]] == false && min_z_triangle[current_index] >= min_z_triangle[save_current_index[i]]) {
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						if (candidate_triangles[current_index][j] == remove_triangles[i][k]) {
							//if(height_of_beam_search != 2)
							remove_triangles.push_back(candidate_triangles[current_index]);
							id_of_furcation_of_blocks.push_back(id_of_furcation_of_blocks[i]);
							all_furcation_of_blocks[id_of_furcation_of_blocks[i]].push_back(candidate_triangles[current_index]);
							//////all_furcation_of_blocks[t].push_back(candidate_triangles[current_index]);

							Eigen::Vector3d face_normal = faceNormal(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							face_normal.normalize();
							Eigen::Vector3d base_normal(0, 0, 1);
							int jud_self_support = (face_normal.dot(base_normal) + sin(max_self_support_slope_angle) >= 0);
							double current_area = triangleArea(all_slicer.positions[remove_triangles[remove_triangles.size() - 1][0]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][1]], all_slicer.positions[remove_triangles[remove_triangles.size() - 1][2]]);
							sum_area += current_area;

							if (jud_self_support == false) {
								sum_not_self_support_area += current_area;
								if (sum_not_self_support_area / sum_area > ratio_of_area_threshold && sum_area >= sum_area_threshold) {
									all_value temp_all_value;
									temp_all_value.value_of_self_support = 0;
									return temp_all_value;
								}
							}

							save_current_index.push_back(current_index);
							jud_triangle_have_been_added[id_triangles[current_index]] = true;
							flag_break = true;
							break;
						}
					}
					if (flag_break == true)
						break;
				}
			}
			if (flag_break == true)
				break;
		}
		current_index++;
		if (current_index >= candidate_triangles.size())
			break;
	}end_time = clock();
	//add remaining face, set a distance threshold Dis, only a face less Dis from other cut layers and no other dependency layer exist, the face is consider to remaining face
	/*double Dis = dh * 4;
	for (int i = 0; i < jud_triangle_have_been_added.size(); i++) {
		if (jud_triangle_have_been_added[i] == false) {
			int id_layer;
			double min_dis = 99999999;
			for (int t = 0; t < all_cut_layers.size(); t++) {
				if (all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] <= 2 * dh && all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] >= -0.001) {
					for (int j = 0; j < all_cut_layers[t].size(); j++) {
						cv::Point3d current_triangle_point(slicer.positions[slicer.triangles[i][0]][0], slicer.positions[slicer.triangles[i][0]][1], slicer.positions[slicer.triangles[i][0]][2]);
						cv::Point3d current_layer_point(all_cut_layers[t][j].x, all_cut_layers[t][j].y, all_cut_layers[t][0].z);

						double distance = distance_3d(current_triangle_point, current_layer_point);
						if (distance < min_dis) {
							min_dis = distance;
							id_layer = t;
						}
					}
				}
			}
			if (min_dis < Dis && all_cut_layers_dependency_layer[id_layer] == 0)
				remove_triangles.push_back(slicer.triangles[i]);
		}
	}*/

	all_value all_calculated_value;
	all_calculated_value.number_of_remaining_face = all_slicer.triangles.size() - remove_triangles.size();
	all_slicer.triangles = remove_triangles;


	all_calculated_value.value_of_projected = calculate_projected_area(all_slicer, all_furcation_of_blocks, all_cut_layers);

	double min_z = 9999999, max_z = -9999999;
	double value_of_manufacturing_require;
	for (int i = 0; i < all_slicer.positions.size(); i++) {
		min_z = min(min_z, all_slicer.positions[i][2]);
		max_z = max(max_z, all_slicer.positions[i][2]);
	}
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		double current_area = triangleArea(all_slicer.positions[all_slicer.triangles[i][0]], all_slicer.positions[all_slicer.triangles[i][1]], all_slicer.positions[all_slicer.triangles[i][2]]);
		value_of_manufacturing_require = vasco::calculateManufacturingRequire(all_slicer, i, min_z, max_z);
		sum_value_of_manufacturing_require += value_of_manufacturing_require * current_area;
	}

	all_calculated_value.value_of_self_support = 1 - sum_not_self_support_area / sum_area;
	all_calculated_value.large_base = sum_value_of_manufacturing_require / sum_area;


	//将slicer旋转回原始位置
	RotateSlicerPositions(all_slicer, vector_after, vectorBefore);

	auto filtered_patch = FilterSurfaceRemoveTriangles(all_slicer, remove_triangles);
	if (!all_slicer.positions.empty() && !filtered_patch.empty()) {
		const std::string mesh_name = "all_slicer_" + std::to_string(height_of_beam_search) + "_" + std::to_string(cont_number_of_queue);
		polyscope::registerSurfaceMesh(mesh_name, all_slicer.positions, filtered_patch);
	}
	//polyscope::show();

	//cout << "()()()(" << double(end_time - start_time) / CLOCKS_PER_SEC << endl;
	return all_calculated_value;
}

void HybridManufacturing::SortCutLayersByHeight(
	vector<vector<cv::Point3d>>& all_cut_layers,
	vector<int>& all_cut_layers_dependency_layer,
	vector<int>& flag_cut_layers_is_hole,
	map<int, int>& follow_index) const
{
	follow_index.clear();
	for (int i = 0; i < all_cut_layers.size(); i++)
		follow_index.insert({ i,i });	//初始化follow_index
	for (int i = 0; i < all_cut_layers.size(); i++) {
		for (int j = i + 1; j < all_cut_layers.size(); j++) {	//按照z值从小到大排序cut layers
			if (all_cut_layers[i][0].z > all_cut_layers[j][0].z) {
				swap(all_cut_layers[i], all_cut_layers[j]);	//交换cut layers
				swap(follow_index[i], follow_index[j]);	//交换对应关系
				swap(flag_cut_layers_is_hole[i], flag_cut_layers_is_hole[j]);	//交换是否为孔洞的标记
				swap(all_cut_layers_dependency_layer[i], all_cut_layers_dependency_layer[j]);	//交换对应的依赖层数量
			}
		}
	}
}

Slicer_2 HybridManufacturing::LoadSlicerForCutMesh(
	bool flag_is_continue_block,
	int height_of_beam_search,
	int index_of_pre_node,
	int id_continue) const
{
	Slicer_2 slicer;
	if (!flag_is_continue_block) {
		slicer.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + ".obj");
		cout << "&" << endl;
	}
	else {
		slicer.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + "_" + to_string(id_continue - 1) + "_subblock.obj");
		cout << "*" << endl;
	}
	return slicer;
}

void HybridManufacturing::RotateSlicerPositions(
	Slicer_2& slicer,
	const Eigen::Vector3d& vector_before,
	const Eigen::Vector3d& vector_after) const
{
	const Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_before, vector_after).toRotationMatrix();
	for (int i = 0; i < slicer.positions.size(); i++) {
		Eigen::MatrixXd temp_V(3, 1);
		temp_V(0, 0) = slicer.positions[i][0];
		temp_V(1, 0) = slicer.positions[i][1];
		temp_V(2, 0) = slicer.positions[i][2];
		temp_V = rotMatrix.inverse() * temp_V;
		slicer.positions[i][0] = temp_V(0, 0);
		slicer.positions[i][1] = temp_V(1, 0);
		slicer.positions[i][2] = temp_V(2, 0);
	}
}

void HybridManufacturing::RotateLayersForVisualization(
	vector<vector<cv::Point3d>>& all_layers,
	vector<vector<cv::Point3d>>& all_layers_contain,
	const Eigen::Vector3d& vector_after,
	const Eigen::Vector3d& vector_before) const
{
	const Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vector_after, vector_before).toRotationMatrix();
	for (int i = 0; i < all_layers.size(); i++) {
		for (int j = 0; j < all_layers[i].size(); j++) {
			Eigen::MatrixXd temp_V(3, 1);
			temp_V(0, 0) = all_layers[i][j].x;
			temp_V(1, 0) = all_layers[i][j].y;
			temp_V(2, 0) = all_layers[i][j].z;
			temp_V = rotMatrix.inverse() * temp_V;
			all_layers[i][j].x = temp_V(0, 0);
			all_layers[i][j].y = temp_V(1, 0);
			all_layers[i][j].z = temp_V(2, 0);
		}

		for (int j = 0; j < all_layers_contain[i].size(); j++) {
			Eigen::MatrixXd temp_V(3, 1);
			temp_V(0, 0) = all_layers_contain[i][j].x;
			temp_V(1, 0) = all_layers_contain[i][j].y;
			temp_V(2, 0) = all_layers_contain[i][j].z;
			temp_V = rotMatrix.inverse() * temp_V;
			all_layers_contain[i][j].x = temp_V(0, 0);
			all_layers_contain[i][j].y = temp_V(1, 0);
			all_layers_contain[i][j].z = temp_V(2, 0);
		}
	}
}

void HybridManufacturing::VisualizeCutLayers(
	const vector<vector<cv::Point3d>>& all_layers,
	const vector<vector<cv::Point3d>>& all_layers_contain,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	bool judge_continue_additive,
	int id_continue,
	const Eigen::Vector3d& vector_after) const
{
	Visual Vis;
	cout << id_continue << endl;
	Vis.generateModelForRendering_5(all_layers, all_layers_contain, height_of_beam_search, cont_number_of_queue, file_name, index_of_pre_node, judge_continue_additive, id_continue);
	string vis_file(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "(" + to_string(index_of_pre_node) + ")_Layer");
	string vis_file_contain(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "(" + to_string(index_of_pre_node) + ")_Layer_contain");
	Eigen::Vector3f current_orientation = vector_after.cast<float>();
	if (open_vis_stair_case == true)
		visualize_layers_stair_case(vis_file, vis_file_contain, current_orientation, judge_continue_additive, id_continue);
}

std::vector<std::vector<int>> HybridManufacturing::EvaluateMergedPatchToolCollision(
	const Slicer_2& merged_patch,
	const std::vector<int>& merged_face_source_patch_id,
	cutter cutting_tool) const
{
	const int face_count = static_cast<int>(merged_patch.triangles.size());
	const int ori_count = static_cast<int>(sampling_subtractive.sample_points.size());

	// face x orientation: -1 表示无碰撞；>=0 表示碰撞面的最大 patch_index
	std::vector<std::vector<int>> min_collision_patch_matrix(face_count, std::vector<int>(ori_count, -1));
	if (face_count == 0 || ori_count == 0) {
		return min_collision_patch_matrix;
	}

	// 每个三角面的采样点（内心）
	std::vector<vasco::core::Vec3> sample_points(face_count);
	for (int i = 0; i < face_count; ++i) {
		const auto& tri = merged_patch.triangles[i];
		const auto& v1 = merged_patch.positions[tri[0]];
		const auto& v2 = merged_patch.positions[tri[1]];
		const auto& v3 = merged_patch.positions[tri[2]];

		double a = distance3d(v1, v2);
		double b = distance3d(v1, v3);
		double c = distance3d(v2, v3);

		vasco::core::Vec3 p{};
		if (a + b + c == 0.0) {
			p = v1;
		}
		else {
			for (int k = 0; k < 3; ++k) {
				p[k] = (a * v1[k] + b * v2[k] + c * v3[k]) / (a + b + c);
			}
		}
		sample_points[i] = p;
	}

	PrepareToolForCollision(cutting_tool);

	for (int ori = 0; ori < ori_count; ++ori) {
		std::vector<std::vector<Eigen::MatrixXd>> temp_faces(face_count, std::vector<Eigen::MatrixXd>(3));
		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				temp_faces[i][k].resize(3, 1);
				int vid = merged_patch.triangles[i][k];
				temp_faces[i][k](0, 0) = merged_patch.positions[vid][0];
				temp_faces[i][k](1, 0) = merged_patch.positions[vid][1];
				temp_faces[i][k](2, 0) = merged_patch.positions[vid][2];
			}
		}

		std::vector<Eigen::MatrixXd> temp_samples(face_count);
		for (int i = 0; i < face_count; ++i) {
			temp_samples[i].resize(3, 1);
			temp_samples[i](0, 0) = sample_points[i][0];
			temp_samples[i](1, 0) = sample_points[i][1];
			temp_samples[i](2, 0) = sample_points[i][2];
		}

		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();

		Eigen::Matrix3d rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();

		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				temp_faces[i][k] = rotMatrix.inverse() * temp_faces[i][k];
			}
			temp_samples[i] = rotMatrix.inverse() * temp_samples[i];
		}

		std::vector<double> max_z_of_faces(face_count, MIN_D);
		for (int i = 0; i < face_count; ++i) {
			for (int k = 0; k < 3; ++k) {
				max_z_of_faces[i] = std::max(max_z_of_faces[i], temp_faces[i][k](2, 0));
			}
		}

		std::vector<Eigen::Vector3d> normals(face_count);
		for (int i = 0; i < face_count; ++i) {
			const Eigen::Vector3d v1(temp_faces[i][0](0, 0), temp_faces[i][0](1, 0), temp_faces[i][0](2, 0));
			const Eigen::Vector3d v2(temp_faces[i][1](0, 0), temp_faces[i][1](1, 0), temp_faces[i][1](2, 0));
			const Eigen::Vector3d v3(temp_faces[i][2](0, 0), temp_faces[i][2](1, 0), temp_faces[i][2](2, 0));

			Eigen::Vector3d n = (v2 - v1).cross(v3 - v1);
			n.normalize();
			normals[i] = n;
		}

		// 计算：刀尖在第 i 面采样点、方向 ori 时，碰撞到的最大 patch_index
		for (int i = 0; i < face_count; ++i) {
			Eigen::Vector3d center_point;
			center_point.x() = temp_samples[i](0, 0) + cutting_tool.cylinder_r * normals[i].x();
			center_point.y() = temp_samples[i](1, 0) + cutting_tool.cylinder_r * normals[i].y();
			center_point.z() = temp_samples[i](2, 0) + cutting_tool.cylinder_r * normals[i].z();

			int max_patch_idx = std::numeric_limits<int>::min();
			bool has_collision = false;

			for (int ii = 0; ii < face_count; ++ii) {
				if (ii == i) {
					continue;
				}
				if (CheckToolCollisionWithCell(center_point, temp_faces[ii], max_z_of_faces[ii], cutting_tool, 30.0, 3.0)) {
					has_collision = true;
					if (ii < static_cast<int>(merged_face_source_patch_id.size())) {
						max_patch_idx = std::max(max_patch_idx, merged_face_source_patch_id[ii]);
					}
				}
			}

			// 无碰撞记为 -1；有碰撞则赋值为patch index + 1
			if (!has_collision) {
				min_collision_patch_matrix[i][ori] = -1;
			}
			else {
				min_collision_patch_matrix[i][ori] = max_patch_idx + 1;
			}
		}
	}

	std::vector<int> min_collision_patch_per_face(face_count, -1);
	for (int i = 0; i < face_count; ++i) {
		int row_min = std::numeric_limits<int>::max();
		for (int ori = 0; ori < ori_count; ++ori) {
			row_min = std::min(row_min, min_collision_patch_matrix[i][ori]);
		}
		if (row_min == std::numeric_limits<int>::max()) {
			row_min = -1;
		}
		min_collision_patch_per_face[i] = row_min;
	}

	// 可选输出
	std::cout << "[Info] min_collision_patch_per_face size="
		<< min_collision_patch_per_face.size() << std::endl;

	// 根据 min_collision_patch_per_face 输出带颜色的 OBJ（每个面一个颜色）
	ExportMergedPatchFaceColorOBJ(
		merged_patch,
		min_collision_patch_per_face,
		".\\vis\\merged_patch_face_color_by_min_collision_patch.obj");

	return min_collision_patch_matrix;
}

void HybridManufacturing::ExportMergedPatchFaceColorOBJ(
	const Slicer_2& merged_patch,
	const std::vector<int>& max_collision_patch_per_face,
	const std::string& color_obj_file) const
{
	std::ofstream ofs(color_obj_file);
	if (!ofs.is_open()) {
		std::cout << "[Warn] cannot open file for writing: " << color_obj_file << std::endl;
		return;
	}

	int max_nonneg_label = -1;
	for (int v : max_collision_patch_per_face) {
		if (v >= 0) {
			max_nonneg_label = std::max(max_nonneg_label, v);
		}
	}
	std::cout << "[EvaluateMergedPatchToolCollision] max_nonneg_label = "
		<< max_nonneg_label << std::endl;

	auto color_from_label = [max_nonneg_label](int label) -> std::array<double, 3> {
		// -1: 无碰撞，固定颜色（浅灰）
		if (label < 0) {
			return { 0.80, 0.80, 0.80 };
		}

		// 非负: 越大越深（蓝色系）
		double t = 0.0;
		if (max_nonneg_label > 0) {
			t = static_cast<double>(label) / static_cast<double>(max_nonneg_label);
		}
		t = std::max(0.0, std::min(1.0, t));

		const std::array<double, 3> light = { 0.78, 0.88, 1.00 };
		const std::array<double, 3> dark = { 0.05, 0.20, 0.60 };

		std::array<double, 3> c;
		for (int k = 0; k < 3; ++k) {
			c[k] = light[k] * (1.0 - t) + dark[k] * t;
		}
		return c;
		};

	// 每个面写3个独立顶点，保证“面颜色”生效
	int v_count = 0;
	for (int i = 0; i < static_cast<int>(merged_patch.triangles.size()); ++i) {
		const auto& tri = merged_patch.triangles[i];
		const auto color = color_from_label(max_collision_patch_per_face[i]);

		for (int k = 0; k < 3; ++k) {
			const auto& p = merged_patch.positions[tri[k]];
			ofs << "v " << p[0] << " " << p[1] << " " << p[2] << " "
				<< color[0] << " " << color[1] << " " << color[2] << "\n";
		}

		ofs << "f " << (v_count + 1) << " " << (v_count + 2) << " " << (v_count + 3) << "\n";
		v_count += 3;
	}

	std::cout << "[Info] wrote colored OBJ: " << color_obj_file << std::endl;
}

Slicer_2 HybridManufacturing::MergeBlockPatchesWithDedup(
	int max_patch_index,
	std::vector<int>& merged_vertex_source_patch_id,
	std::vector<int>& merged_face_source_patch_id,
	double merge_eps) const
{
	Slicer_2 merged;
	merged_vertex_source_patch_id.clear();
	merged_face_source_patch_id.clear();

	struct QuantizedVertexKey {
		long long x;
		long long y;
		long long z;
		bool operator<(const QuantizedVertexKey& rhs) const {
			if (x != rhs.x) return x < rhs.x;
			if (y != rhs.y) return y < rhs.y;
			return z < rhs.z;
		}
	};

	const double inv_eps = 1.0 / merge_eps;
	auto make_vertex_key = [inv_eps](const vasco::Slicer::Vec3& p) {
		return QuantizedVertexKey{
			static_cast<long long>(std::llround(p[0] * inv_eps)),
			static_cast<long long>(std::llround(p[1] * inv_eps)),
			static_cast<long long>(std::llround(p[2] * inv_eps))
		};
		};

	std::map<QuantizedVertexKey, int> vertex_map;
	std::map<std::array<int, 3>, int> face_key_to_merged_face_index;

	for (int patch_index = 1; patch_index <= max_patch_index; ++patch_index) {
		Slicer_2 patch;
		const std::string patch_file = ".\\vis\\block_patch-" + to_string(patch_index) + "_" + ".obj";

		if (!patch.load(patch_file)) {
			std::cout << "[Warn] cannot load patch file: " << patch_file << std::endl;
			continue;
		}

		std::vector<int> local_to_global(patch.positions.size(), -1);

		for (int i = 0; i < static_cast<int>(patch.positions.size()); ++i) {
			const auto key = make_vertex_key(patch.positions[i]);
			auto it = vertex_map.find(key);

			if (it == vertex_map.end()) {
				const int new_index = static_cast<int>(merged.positions.size());
				merged.positions.push_back(patch.positions[i]);
				vertex_map.insert({ key, new_index });
				local_to_global[i] = new_index;

				// 单归属：首次出现即最小 patch_index
				merged_vertex_source_patch_id.push_back(patch_index);
			}
			else {
				// 已存在：保持原归属（更小 patch_index）
				local_to_global[i] = it->second;
			}
		}

		for (const auto& tri : patch.triangles) {
			int a = local_to_global[tri[0]];
			int b = local_to_global[tri[1]];
			int c = local_to_global[tri[2]];

			if (a == b || b == c || a == c) {
				continue;
			}

			std::array<int, 3> key = { a, b, c };
			std::sort(key.begin(), key.end());

			auto fit = face_key_to_merged_face_index.find(key);
			if (fit == face_key_to_merged_face_index.end()) {
				const int new_face_index = static_cast<int>(merged.triangles.size());
				merged.triangles.push_back({ a, b, c });
				face_key_to_merged_face_index.insert({ key, new_face_index });

				// 单归属：首次出现即最小 patch_index
				merged_face_source_patch_id.push_back(patch_index);
			}
			// 已存在：保持原归属（更小 patch_index）
		}
	}

	return merged;
}


void HybridManufacturing::CutMesh(
	CutLayerVector all_layers,
	CutLayerVector all_layers_contain,
	CutLayerVector all_cut_layers,
	Eigen::Vector3d vector_after,
	int height_of_beam_search,
	int cont_number_of_queue,
	int index_of_pre_node,
	vector<int> all_cut_layers_dependency_layer,
	bool& jud_outer_beam_search_terminate,
	vector<TRiangle>& current_remove_triangles,
	Slicer_2& current_slicer,
	bool judge_continue_additive,
	bool flag_is_continue_block,
	int pre_cont_number_of_queue,
	vector<bool>& jud_error,
	int id_node,
	int id_continue,
	vector<int> flag_cut_layers_is_hole)
{
	bool using_solid_model = true;

	std::vector<bool> jud_triangle_have_been_added;	//三角形是否被添加进remove_triangles
	vector<vector<int>> cutting_plane_points;	//每个切割平面对应的切割点
	cutting_plane_points.resize(all_cut_layers.size());
	vector<vector<pair<int, int>>> cutting_plane_edges;	//每个切割平面对应的切割边 <顶点index,顶点index>
	cutting_plane_edges.resize(all_cut_layers.size());

	if (all_cut_layers.size() == 0) {
		std::cout << "[HybridManufacturing::CutMesh]No cut layers!" << std::endl;
		return;
	}
	//////sort cut layers//////
	map<int, int> follow_index;	//记录排序前后cut layers的对应关系 <排序后index,排序前index>
	SortCutLayersByHeight(all_cut_layers, all_cut_layers_dependency_layer, flag_cut_layers_is_hole, follow_index);
	///////////////////////////
	//cout << "a" << id_continue;
	Slicer_2 slicer = LoadSlicerForCutMesh(flag_is_continue_block, height_of_beam_search, index_of_pre_node, id_continue);


	if (flag_is_continue_block == true) {
		height_of_beam_search--;
		cont_number_of_queue = pre_cont_number_of_queue;
	}

	//need rotate first//
	Eigen::Vector3d vectorBefore(0, 0, 1);
	RotateSlicerPositions(slicer, vectorBefore, vector_after);
	/////////////////////

	//layer visualization//
	RotateLayersForVisualization(all_layers, all_layers_contain, vector_after, vectorBefore);
	cout << "b" << endl;
	VisualizeCutLayers(all_layers, all_layers_contain, height_of_beam_search, cont_number_of_queue, index_of_pre_node, judge_continue_additive, id_continue, vector_after);
	////////////////////////
	cout << "B" << endl;
	//--------------------cut-----------------------//
	slicer.normal[0] = 0;
	slicer.normal[1] = 0;
	slicer.normal[2] = 1;

	slicer.jud_plane.resize(slicer.triangles.size());
	for (int i = 0; i < slicer.triangles.size(); i++)
		slicer.jud_plane[i] = false;
	clock_t start_time_3, end_time_3;

	for (int i = 0; i < all_cut_layers.size(); i++) {
		slicer.origin[0] = 0;
		slicer.origin[1] = 0;
		slicer.origin[2] = all_cut_layers[i][0].z;
		slicer.cut();
		//cout << "ok";
		cout << "i " << i << " slicer.origin[2] " << slicer.origin[2] << endl;
	}

	std::vector<TRiangle> ori_triangle = slicer.triangles;	//保存切割后的所有三角形
	jud_triangle_have_been_added.resize(slicer.triangles.size());	//为ori_triangle中的每个三角形分配一个标记，表示是否被添加进remove_triangles
	for (int i = 0; i < slicer.triangles.size(); i++)
		jud_triangle_have_been_added[i] = false;

	Slicer_2 all_slicer;	//all_slicer保存切割后的所有三角形和顶点

	all_slicer.positions = slicer.positions;
	all_slicer.triangles = ori_triangle;
	//--------------------save candidate_triangles-----------------------//
	int current_index = 0;	//
	std::vector<VEctor> min_z_point;	//每个候选三角形中z值最小的顶点坐标
	std::vector<double> min_z_triangle;	//每个候选三角形中z值最小的顶点的z值
	vector<int> index_of_min_point_in_triangle;	//每个候选三角形中z值最小的顶点在三角形中的索引
	std::vector<TRiangle> candidate_triangles;	//候选三角形集合
	std::vector<int> id_candidate_triangles;	//候选三角形在all_slicer.triangles中的索引
	std::vector<int> id_triangles;	//候选三角形在slicer.triangles中的索引
	std::vector<TRiangle> remove_triangles;	//最终被移除的三角形集合
	std::vector<int> id_remove_triangles;	//最终被移除的三角形在all_slicer.triangles中的索引
	candidate_triangles.clear();
	id_triangles.clear();
	double min_z_all_cut_layers = 9999999;	//所有切割层中z值最小的值
	for (int t = 0; t < all_cut_layers.size(); t++)
		min_z_all_cut_layers = min(min_z_all_cut_layers, all_cut_layers[t][0].z);
	start_time_3 = clock();
	for (int i = 0; i < all_slicer.triangles.size(); i++) {	//遍历所有切割后的三角形
		double temp_min_z_triangle = 9999999;	//当前三角形中z值最小的顶点的z值
		for (int k = 0; k < 3; k++) {
			temp_min_z_triangle = min(all_slicer.positions[all_slicer.triangles[i][k]][2], temp_min_z_triangle);
		}
		if (temp_min_z_triangle + 0.001 >= min_z_all_cut_layers) {	//如果当前三角形中z值最小的顶点的z值大于等于所有切割层中z值最小的值，则将该三角形加入候选三角形集合
			candidate_triangles.push_back(all_slicer.triangles[i]);	//加入候选三角形集合
			id_candidate_triangles.push_back(i);	//记录该候选三角形在all_slicer.triangles中的索引
			for (int k = i; k < slicer.triangles.size(); k++)
				if (slicer.triangles[k] == all_slicer.triangles[i]) {	//找到该候选三角形在slicer.triangles中的索引
					id_triangles.push_back(k);
					break;
				}
		}
	}
	end_time_3 = clock();
	for (int i = 0; i < candidate_triangles.size(); i++) {
		double temp_min_z_triangle = 9999999;	//当前候选三角形中z值最小的顶点的z值
		VEctor temp_min_z_point = all_slicer.positions[candidate_triangles[i][0]];
		int temp_index_of_min_point_in_triangle;
		for (int k = 0; k < 3; k++) {
			if (all_slicer.positions[candidate_triangles[i][k]][2] < temp_min_z_triangle) {	//找到当前候选三角形中z值最小的顶点
				temp_min_z_triangle = all_slicer.positions[candidate_triangles[i][k]][2];	//更新z值最小的顶点的z值
				temp_min_z_point = all_slicer.positions[candidate_triangles[i][k]];	//更新z值最小的顶点的坐标
				temp_index_of_min_point_in_triangle = k;	//更新z值最小的顶点在三角形中的索引
			}
		}
		min_z_triangle.push_back(temp_min_z_triangle);	//将当前候选三角形中z值最小的顶点的z值加入min_z_triangle
		min_z_point.push_back(temp_min_z_point);	//将当前候选三角形中z值最小的顶点的坐标加入min_z_point
		index_of_min_point_in_triangle.push_back(temp_index_of_min_point_in_triangle);	//将当前候选三角形中z值最小的顶点在三角形中的索引加入index_of_min_point_in_triangle
	}

	cout << "()()()(" << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << endl;
	clock_t start_time_2, end_time_2;
	start_time_2 = clock();
	for (int i = 0; i < candidate_triangles.size(); i++) {	//对候选三角形按照min_z_triangle进行排序
		for (int j = i + 1; j < candidate_triangles.size(); j++) {
			if (min_z_triangle[i] > min_z_triangle[j]) {
				swap(candidate_triangles[i], candidate_triangles[j]);	//交换候选三角形
				swap(id_candidate_triangles[i], id_candidate_triangles[j]);	//交换候选三角形在all_slicer.triangles中的索引
				swap(id_triangles[i], id_triangles[j]);	//交换候选三角形在slicer.triangles中的索引
				swap(min_z_triangle[i], min_z_triangle[j]);	//交换min_z_triangle
				swap(min_z_point[i], min_z_point[j]);	//交换min_z_point
				swap(index_of_min_point_in_triangle[i], index_of_min_point_in_triangle[j]);	//交换index_of_min_point_in_triangle
			}
		}
	}
	end_time_2 = clock();
	cout << "()()()(" << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << endl;
	//slicer.triangles = candidate_triangles;

	//--------------------save OPP_triangles one by one-----------------------//
	vector<int> save_current_index;
	for (int t = 0; t < all_cut_layers.size(); t++) {	//枚举all_cut_layers中的每一层
		current_index = 0;   //该方式也许比较慢
		double boundary_bottom = 999999, boundary_left = 999999, boundary_top = -999999, boundary_right = -999999;	//当前切割层的边界
		//Point_2* points = new Point_2[all_cut_layers[t].size()];
		vector<Point_2> points;
		points.resize(all_cut_layers[t].size());
		for (int i = 0; i < all_cut_layers[t].size(); i++) {	//枚举all_cut_layers[t]中的点，计算当前切割层的二维AABB
			boundary_top = std::max(boundary_top, all_cut_layers[t][i].y);
			boundary_bottom = std::min(boundary_bottom, all_cut_layers[t][i].y);
			boundary_right = std::max(boundary_right, all_cut_layers[t][i].x);
			boundary_left = std::min(boundary_left, all_cut_layers[t][i].x);
			Point_2 temp_point(all_cut_layers[t][i].x, all_cut_layers[t][i].y);
			points[i] = temp_point;
		}
		cv::Point2d current_triangle_point;
		cv::Point2d current_layer_point;
		for (; current_index < candidate_triangles.size(); current_index++) {	//枚举所有候选三角形
			if (abs(min_z_triangle[current_index] - all_cut_layers[t][0].z) > 0.0001) {
				if (min_z_triangle[current_index] > all_cut_layers[t][0].z) {	//如果当前候选三角形中z值最小的顶点的z值大于当前切割层的z值，则跳出循环
					break;
				}
			}
			else {	//如果当前候选三角形中最小z值等于当前切割层的z值,才进行是否加入remove_triangles的判断
				bool jud_is_boundary_point = false;
				//for (int i = 0; i < all_cut_layers[t].size(); i++) {
				int cont_inside_boundary = 0;
				for (int k = 0; k < 3; k++) {	//枚举当前候选三角形的三个顶点，判断是否在当前切割层的AABB内
					if (all_slicer.positions[candidate_triangles[current_index][k]][0] + 0.01 >= boundary_left && all_slicer.positions[candidate_triangles[current_index][k]][0] - 0.01 <= boundary_right
						&& all_slicer.positions[candidate_triangles[current_index][k]][1] + 0.01 >= boundary_bottom && all_slicer.positions[candidate_triangles[current_index][k]][1] - 0.01 <= boundary_top) {
						cont_inside_boundary++;	//判断当前顶点在AABB内，计数加1
					}
				}
				if (cont_inside_boundary >= 2)
					jud_is_boundary_point = true;
				/*if (min_z_point[current_index][0] + 0.01 >= boundary_left && min_z_point[current_index][0] - 0.01 <= boundary_right &&
					min_z_point[current_index][1] + 0.01 >= boundary_bottom && min_z_point[current_index][1] - 0.01 <= boundary_top) {
					jud_is_boundary_point = true;
					break;
				}*/
				//}

				if (jud_is_boundary_point == true) {
					auto layer_begin = points.begin();
					auto layer_end = points.end();
					//if (check_inside_2(Point_2(current_triangle_point.x, current_triangle_point.y), points, points + all_cut_layers[t].size(), K())) {
					current_triangle_point.x = min_z_point[current_index][0];
					current_triangle_point.y = min_z_point[current_index][1];
					/*if (check_inside_2(Point_2(current_triangle_point.x, current_triangle_point.y), layer_begin, layer_end, K())) {*/

					if (1) {
						int cont_cutting_points = 0;
						int temp_left = 0, temp_right = 0;
						for (int k = 0; k < 3; k++) {
							if (abs(all_cut_layers[t][0].z - all_slicer.positions[candidate_triangles[current_index][k]][2]) < 0.0001) {	//判断当前顶点是否在切割平面上（前面已经保证了z值最小的顶点不大于切割平面z值）
								cont_cutting_points++;
								if (cont_cutting_points == 1)
									temp_left = candidate_triangles[current_index][k];
								else if (cont_cutting_points == 2) {
									temp_right = candidate_triangles[current_index][k];
									cutting_plane_edges[t].push_back(make_pair(temp_left, temp_right));
									break;
								}
								//cutting_plane_points[t].push_back(candidate_triangles[current_index][index_of_min_point_in_triangle[current_index]]);
							}
						}
						remove_triangles.push_back(candidate_triangles[current_index]);	//将当前候选三角形加入remove_triangles
						id_remove_triangles.push_back(id_candidate_triangles[current_index]);	//将当前候选三角形在all_slicer.triangles中的索引加入id_remove_triangles
						/*if (height_of_beam_search == 2 && id_triangles[current_index] == 10297)
							cout << "********************************** " << all_cut_layers[t][0].z << endl;*/

						save_current_index.push_back(current_index);	//记录当前候选三角形在candidate_triangles中的索引到save_current_index
						jud_triangle_have_been_added[id_triangles[current_index]] = true;	//jud_triangle_have_been_added标记当前候选三角形已被添加进remove_triangles
					}
					else {
						for (int j = 0; j < all_cut_layers[t].size(); j++) {
							current_layer_point.x = all_cut_layers[t][j].x;	//
							current_layer_point.y = all_cut_layers[t][j].y;
							current_triangle_point.x = min_z_point[current_index][0];
							current_triangle_point.y = min_z_point[current_index][1];
							//if ((jud_is_boundary_point == false && distance_2d(current_layer_point, current_triangle_point) < 0.002) || (jud_is_boundary_point == true)) {   //0.002
							if (distance2d(current_layer_point, current_triangle_point) < 0.1) {  //4.0
								int cont_cutting_points = 0;
								int temp_left = 0, temp_right = 0;
								for (int k = 0; k < 3; k++) {
									if (abs(all_cut_layers[t][0].z - all_slicer.positions[candidate_triangles[current_index][k]][2]) < 0.0001) {
										cont_cutting_points++;
										if (cont_cutting_points == 1)
											temp_left = candidate_triangles[current_index][k];
										else if (cont_cutting_points == 2) {
											temp_right = candidate_triangles[current_index][k];
											cutting_plane_edges[t].push_back(make_pair(temp_left, temp_right));
											break;
										}
										//cutting_plane_points[t].push_back(candidate_triangles[current_index][index_of_min_point_in_triangle[current_index]]);
									}
								}
								remove_triangles.push_back(candidate_triangles[current_index]);
								id_remove_triangles.push_back(id_candidate_triangles[current_index]);
								/*if (height_of_beam_search == 2 && id_triangles[current_index] == 10297)
									cout << "********************************** " << all_cut_layers[t][0].z << endl;*/

								save_current_index.push_back(current_index);
								jud_triangle_have_been_added[id_triangles[current_index]] = true;
								break;
							}
						}
					}
				}
			}

		}
	}
	current_index = 0;

	while (1) {
		bool flag_break = false;
		for (int i = 0; i < remove_triangles.size(); i++)
		{
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					if (jud_triangle_have_been_added[id_triangles[current_index]] == false && candidate_triangles[current_index][j] == remove_triangles[i][k]
						&& min_z_triangle[current_index] >= min_z_triangle[save_current_index[i]]) {	//如果当前候选三角形的某个顶点与remove_triangles[i]的某个顶点相同，且当前候选三角形中z值最小的顶点的z值大于等于remove_triangles[i]中z值最小的顶点的z值
						//if (height_of_beam_search != 2)
						remove_triangles.push_back(candidate_triangles[current_index]);	//将当前候选三角形加入remove_triangles
						id_remove_triangles.push_back(id_candidate_triangles[current_index]);
						save_current_index.push_back(current_index);
						jud_triangle_have_been_added[id_triangles[current_index]] = true;
						flag_break = true;
						break;
					}
				}
				if (flag_break == true)
					break;
			}
			if (flag_break == true)
				break;
		}
		current_index++;
		if (current_index >= candidate_triangles.size())
			break;
	}

	vasco::Slicer temp_slicer_1;
	temp_slicer_1.clear();
	temp_slicer_1.positions = all_slicer.positions;
	for (int i = 0; i < id_remove_triangles.size(); i++) {
		temp_slicer_1.triangles.push_back(all_slicer.triangles[id_remove_triangles[i]]);
	}
	//cout << "remove faces " << temp_slicer_1.triangles.size() << endl;
	temp_slicer_1.save(file_name + "-all2-" + std::to_string(height_of_beam_search) + "QWERQ.obj");

	//add remaining face, set a distance threshold Dis, only a face less Dis from other cut layers and no other dependency layer exist, the face is consider to remaining face
	double Dis = dh * 2;
	for (int i = 0; i < jud_triangle_have_been_added.size(); i++) {
		if (jud_triangle_have_been_added[i] == false) {
			int id_layer;
			double min_dis = 99999999;
			for (int t = 0; t < all_cut_layers.size(); t++) {
				if (all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] <= 2 * dh && all_cut_layers[t][0].z - slicer.positions[slicer.triangles[i][0]][2] >= -0.001) {	//如果当前切割层的z值与当前三角形中z值最小的顶点的z值之差在[-0.001,2*dh]范围内
					for (int j = 0; j < all_cut_layers[t].size(); j++) {
						cv::Point3d current_triangle_point(slicer.positions[slicer.triangles[i][0]][0], slicer.positions[slicer.triangles[i][0]][1], slicer.positions[slicer.triangles[i][0]][2]);
						cv::Point3d current_layer_point(all_cut_layers[t][j].x, all_cut_layers[t][j].y, all_cut_layers[t][0].z);

						double distance = distance3d(current_triangle_point, current_layer_point);
						if (distance < min_dis) {
							min_dis = distance;
							id_layer = t;
						}
					}
				}
			}
			if (min_dis < Dis && all_cut_layers_dependency_layer[id_layer] == 0) {
				cout << "还真的能到这个地方啊 - 残余面片 min_dis " << min_dis << "Dis " << Dis << std::endl;
				remove_triangles.push_back(slicer.triangles[i]);
				id_remove_triangles.push_back(i);
				//jud_triangle_have_been_added[i] = true;
			}

		}
	}




	/////////////////////////删除残余面片////////////////////////////
	//建立面片邻接关系
	vector<vector<int>> adjacent_faces(slicer.triangles.size());
	//建立不在remove_triangles中的面片之间的邻接关系
	for (int i = 0; i < slicer.triangles.size(); i++) {
		if (jud_triangle_have_been_added[i] == false) {
			for (int j = i + 1; j < slicer.triangles.size(); j++) {
				if (jud_triangle_have_been_added[j] == false) {
					int cont_same_point = 0;
					for (int k = 0; k < 3; k++) {
						for (int l = 0; l < 3; l++) {
							if (slicer.triangles[i][k] == slicer.triangles[j][l]) {
								cont_same_point++;
								break;
							}
						}
					}
					if (cont_same_point >= 2) {	//如果两个面片有两个以上的公共顶点，则认为它们是邻接的
						adjacent_faces[i].push_back(j);
						adjacent_faces[j].push_back(i);
					}
				}
			}
		}
	}
	//分区
	vector<bool> visited(slicer.triangles.size(), false);
	vector<vector<int>> save_faces(0);
	for (int i = 0; i < slicer.triangles.size(); ++i)
	{
		if (jud_triangle_have_been_added[i] || visited[i]) {
			continue; // 已标记或已访问，跳过
		}

		visited[i] = true;

		// 当前连通分量的索引集合
		std::vector<int> component;
		component.push_back(i);

		// 标准队列进行 BFS
		std::queue<int> q;
		q.push(i);
		while (!q.empty())
		{
			const int u = q.front();
			q.pop();

			// 遍历 u 的邻接面
			const auto& adj = adjacent_faces[u];
			for (int v : adj)
			{
				if (!visited[v] && !jud_triangle_have_been_added[v])
				{
					visited[v] = true;
					q.push(v);
					component.push_back(v);
				}
			}
		}

		// 将本连通分量保存（对应之前的 save_faces.push_back(...)）
		save_faces.push_back(std::move(component));
	}

	// 移除过小的残余面片分量
	for (const auto& comp : save_faces)
	{
		cout << "&&*&& " << comp.size() << endl;
		if (comp.size() < 80)
		{
			for (int idx : comp)
			{
				remove_triangles.push_back(slicer.triangles[idx]);
				id_remove_triangles.push_back(idx);
				jud_triangle_have_been_added[idx] = true;
			}
		}
	}
	//////////////////////////////////////////////////////////////////



	cout << "c" << endl;
	//current_remove_triangles = remove_triangles;

	current_remove_triangles = remove_triangles; // current_remove_triangles先不进行筛选，直接等于remove_triangles，后面再剔除非表面的面片
	cout << "ccc" << endl;

	//将slicer旋转回原始位置
	RotateSlicerPositions(slicer, vector_after, vectorBefore);
	current_slicer = slicer;

	//从all_slicer中移除remove_triangles
	//cout << current_remove_triangles.size() << endl;
	//if (height_of_beam_search != 2 || id_continue != 1)
	for (int i = 0; i < remove_triangles.size(); i++) {
		for (int j = 0; j < all_slicer.triangles.size(); ) {
			if (remove_triangles[i][0] == all_slicer.triangles[j][0] && remove_triangles[i][1] == all_slicer.triangles[j][1] && remove_triangles[i][2] == all_slicer.triangles[j][2]) {
				all_slicer.triangles.erase(all_slicer.triangles.begin() + j);
				break;
			}
			j++;
		}
	}	//潜在问题？：如果顶点顺序不同（例如顺时针/逆时针），将无法删除；可改为排序比较或集合比较。
	//这个潜在问题应该不会出现，因为计算remove_triangles时，是直接从all_slicer.triangles中取出的三角形，三角形顶点顺序没有改变。


	size_t face_cnt = Normals.rows();
	std::cout << "face_cnt " << face_cnt << std::endl;
	current_remove_triangles = FilterSurfaceRemoveTriangles(slicer, remove_triangles);

	//add cutting plane triangles
	/*all_slicer.triangles.insert(all_slicer.triangles.begin(),cutting_plane_points)
	cutting_plane_points*/

	//erase some cutting plane
	for (int i = 0; i < all_cut_layers.size(); i++) {
		if (all_cut_layers_dependency_layer[i] == 0) {
			all_cut_layers.erase(all_cut_layers.begin() + i);
			all_cut_layers_dependency_layer.erase(all_cut_layers_dependency_layer.begin() + i);
			cutting_plane_edges.erase(cutting_plane_edges.begin() + i);
			cutting_plane_points.erase(cutting_plane_points.begin() + i);
			i--;
		}
	}
	cout << "cccc" << endl;
	//sort cutting_plane_points by adjacency relation
	vector<vector<int>> real_cutting_plane_triangles;	//为每个切割层存储排序后的切割平面顶点索引
	real_cutting_plane_triangles.resize(all_cut_layers.size());	//real_cutting_plane_triangles大小与all_cut_layers相同
	for (int i = 0; i < all_cut_layers.size(); i++) {
		cv::Point2d current_triangle_point;
		cv::Point2d current_layer_point;
		current_layer_point.x = all_cut_layers[i][0].x;
		current_layer_point.y = all_cut_layers[i][0].y;
		double min_dis = MAX_D;
		int index_start_point_id;
		for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
			current_triangle_point.x = all_slicer.positions[cutting_plane_edges[i][j].second][0];
			current_triangle_point.y = all_slicer.positions[cutting_plane_edges[i][j].second][1];
			if (distance2d(current_layer_point, current_triangle_point) < min_dis)
			{
				index_start_point_id = j;
				min_dis = distance2d(current_layer_point, current_triangle_point);
			}
		}
		int index_start_point = cutting_plane_edges[i][index_start_point_id].second;
		int index_of_last_edge = index_start_point_id;
		real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][index_start_point_id].first);
		//cout << "kkkk" << cutting_plane_edges[i].size() << endl;
		int jud_select_id = -1;
		int cont_segment = 0;
		while (index_start_point != cutting_plane_edges[i][index_start_point_id].first) {
			bool jud_select = false;
			cont_segment++;
			for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
				if (j != index_of_last_edge && cutting_plane_edges[i][j].first == index_start_point) {
					real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].first);
					index_start_point = cutting_plane_edges[i][j].second;
					index_of_last_edge = j;
					jud_select = true;
					jud_select_id = 0;
					break;
				}
				else if (j != index_of_last_edge && cutting_plane_edges[i][j].second == index_start_point) {
					real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].second);
					index_start_point = cutting_plane_edges[i][j].first;
					index_of_last_edge = j;
					jud_select = true;
					jud_select_id = 1;
					break;
				}
			}
			if (cont_segment > 1000000) {
				jud_select = false;
				//real_cutting_plane_triangles[i].clear();
			}
			if (jud_select == false) {
				//cout << "no:" << i << "  ";
				jud_error[id_node] = true;
				/*for (int j = 0; j < cutting_plane_edges[i].size(); j++) {
					if (j != index_of_last_edge) {
						cv::Point2d current_point_2, current_point_3, current_point_4;
						current_point_3.x = all_slicer.positions[cutting_plane_edges[i][j].first][0];
						current_point_3.y = all_slicer.positions[cutting_plane_edges[i][j].first][1];
						current_point_4.x = all_slicer.positions[cutting_plane_edges[i][j].second][0];
						current_point_4.y = all_slicer.positions[cutting_plane_edges[i][j].second][1];
						if (jud_select_id == 0) {
							current_point_2.x = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].second][0];
							current_point_2.y = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].second][1];
						}
						else if (jud_select_id == 1) {
							current_point_2.x = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].first][0];
							current_point_2.y = all_slicer.positions[cutting_plane_edges[i][index_of_last_edge].first][1];
						}
						if (distance_2d(current_point_2, current_point_3) < 0.00001) {
							real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].first);
							index_start_point = cutting_plane_edges[i][j].second;
							index_of_last_edge = j;
							jud_select = true;
							jud_select_id = 0;
							break;
						}
						else if (distance_2d(current_point_2, current_point_4) < 0.00001) {
							real_cutting_plane_triangles[i].push_back(cutting_plane_edges[i][j].second);
							index_start_point = cutting_plane_edges[i][j].first;
							index_of_last_edge = j;
							jud_select = true;
							jud_select_id = 1;
							break;
						}
					}
				}*/
				break;
			}
		}
	}
	/*if (jud_error[id_node] == true)
		return;*/
		/*cout << "mm" << endl;*/
		/*if (height_of_beam_search == 6)
			return;*/
	for (int i = 0; i < real_cutting_plane_triangles.size(); i++) {
		if (real_cutting_plane_triangles[i].size() == 1) {
			real_cutting_plane_triangles.erase(real_cutting_plane_triangles.begin() + i);
			i--;
		}
	}

	///////////////////////////terminate//////////////////////////////
	if (all_slicer.triangles.size() < terminate_threshold_of_number_of_faces) {
		jud_outer_beam_search_terminate = true;
	}
	//////////////////////////////////////////////////////////////////

	vector<int> id_contact_faces;
	//if (height_of_beam_search != 2 || id_continue != 1)
	if (using_solid_model == true) {
		Anticlockwise(real_cutting_plane_triangles, all_slicer);
		using Coord = double;
		using NN = uint32_t;
		using Point = std::array<Coord, 2>;
		for (int i = 0; i < real_cutting_plane_triangles.size(); i++) {
			if (flag_cut_layers_is_hole[i] != -1)
				continue;
			std::vector<std::vector<Point>> polygon;
			polygon.clear();
			std::vector<Point> temp_vec;
			polygon.push_back(temp_vec);
			map<int, int> map_index_faces;
			map_index_faces.clear();
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				polygon[0].push_back({ all_slicer.positions[real_cutting_plane_triangles[i][j]][0], all_slicer.positions[real_cutting_plane_triangles[i][j]][1] });
				map_index_faces.insert({ j,real_cutting_plane_triangles[i][j] });
			}
			for (int m = 0; m < real_cutting_plane_triangles.size(); m++) {
				if (flag_cut_layers_is_hole[m] != -1 && follow_index[flag_cut_layers_is_hole[m]] == i) {
					for (int j = 0; j < real_cutting_plane_triangles[m].size(); j++) {
						polygon.push_back(temp_vec);
						polygon[1].push_back({ all_slicer.positions[real_cutting_plane_triangles[m][j]][0], all_slicer.positions[real_cutting_plane_triangles[m][j]][1] });
						map_index_faces.insert({ polygon[0].size() + j,real_cutting_plane_triangles[m][j] });
					}
					//多个contour?
					break;
				}
			}

			std::vector<NN> indices = mapbox::earcut<NN>(polygon);

			for (int j = 0; j < indices.size();) {
				TRiangle the_new_cutting_plane_triangle;
				the_new_cutting_plane_triangle[0] = map_index_faces[indices[j]]; j++;
				the_new_cutting_plane_triangle[1] = map_index_faces[indices[j]]; j++;
				the_new_cutting_plane_triangle[2] = map_index_faces[indices[j]]; j++;
				all_slicer.triangles.insert(all_slicer.triangles.end(), the_new_cutting_plane_triangle);  //添加接触面
				id_contact_faces.push_back(all_slicer.triangles.size() - 1);
			}

		}
		//}
	}




	/*if (height_of_beam_search == 1) {
		for (int i = 0; i < real_cutting_plane_triangles.size(); i++) {
			VEctor new_point;
			double temp_x = 0, temp_y = 0;
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				temp_x += all_slicer.positions[real_cutting_plane_triangles[i][j]][0];
				temp_y += all_slicer.positions[real_cutting_plane_triangles[i][j]][1];
			}
			new_point[0] = temp_x / real_cutting_plane_triangles[i].size();
			new_point[1] = temp_y / real_cutting_plane_triangles[i].size();
			new_point[2] = all_slicer.positions[real_cutting_plane_triangles[i][0]][2];
			all_slicer.positions.insert(all_slicer.positions.end(), new_point);
			for (int j = 0; j < real_cutting_plane_triangles[i].size(); j++) {
				TRiangle the_new_cutting_plane_triangle;
				the_new_cutting_plane_triangle[0] = all_slicer.positions.size() - 1;
				the_new_cutting_plane_triangle[1] = real_cutting_plane_triangles[i][j];
				the_new_cutting_plane_triangle[2] = real_cutting_plane_triangles[i][(j + 1) % real_cutting_plane_triangles[i].size()];
				all_slicer.triangles.insert(all_slicer.triangles.end(), the_new_cutting_plane_triangle);
			}
		}
	}*/


	RotateSlicerPositions(all_slicer, vector_after, vectorBefore);

	Slicer_2 all_slicer_2 = all_slicer;
	all_slicer_2.triangles = remove_triangles;
	if (judge_continue_additive == false)
		all_slicer_2.save(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_current" + ".obj");
	else
		all_slicer_2.save(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_current" + "_subblock.obj");

	Eigen::MatrixXd temp_V(all_slicer.positions.size(), 3);
	Eigen::MatrixXi temp_F(all_slicer.triangles.size(), 3);
	for (int i = 0; i < all_slicer.positions.size(); i++) {
		temp_V.row(i).x() = all_slicer.positions[i][0];
		temp_V.row(i).y() = all_slicer.positions[i][1];
		temp_V.row(i).z() = all_slicer.positions[i][2];
	}
	for (int i = 0; i < all_slicer.triangles.size(); i++) {
		temp_F.row(i).x() = all_slicer.triangles[i][0];
		temp_F.row(i).y() = all_slicer.triangles[i][1];
		temp_F.row(i).z() = all_slicer.triangles[i][2];
	}
	if (judge_continue_additive == false) {
		Katana::Instance().stl.saveStlFromObj(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl", temp_V, temp_F);
		vector<std::array<double, 3>> V3;
		vector<std::array<int, 3>> F3;
		vector<std::array<double, 3>> N3;
		ifstream ifs(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl");
		igl::read_stl_ascii(ifs, V3, F3, N3);
		Eigen::MatrixXd V4(V3.size(), 3);
		Eigen::MatrixXi F4(F3.size(), 3);
		Eigen::MatrixXd N4(N3.size(), 3);
		for (int i = 0; i < V3.size(); i++)
			for (int j = 0; j < 3; j++)
				V4(i, j) = V3[i][j];
		for (int i = 0; i < F3.size(); i++)
			for (int j = 0; j < 3; j++)
				F4(i, j) = F3[i][j];
		igl::writeSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_B.stl", V4, F4, igl::FileEncoding::Binary);

		Geometry tessel;
		tessel.visit(ImportSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_B.stl"));
		tessel.visit(ExportOBJ(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".obj"));
	}
	else {
		Katana::Instance().stl.saveStlFromObj(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock.stl", temp_V, temp_F);
		vector<std::array<double, 3>> V3;
		vector<std::array<int, 3>> F3;
		vector<std::array<double, 3>> N3;
		ifstream ifs(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock.stl");
		igl::read_stl_ascii(ifs, V3, F3, N3);
		Eigen::MatrixXd V4(V3.size(), 3);
		Eigen::MatrixXi F4(F3.size(), 3);
		Eigen::MatrixXd N4(N3.size(), 3);
		for (int i = 0; i < V3.size(); i++)
			for (int j = 0; j < 3; j++)
				V4(i, j) = V3[i][j];
		for (int i = 0; i < F3.size(); i++)
			for (int j = 0; j < 3; j++)
				F4(i, j) = F3[i][j];
		igl::writeSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock_B.stl", V4, F4, igl::FileEncoding::Binary);

		Geometry tessel;
		tessel.visit(ImportSTL(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock_B.stl"));
		tessel.visit(ExportOBJ(file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock.obj"));
	}

	string str_contact_faces;
	if (judge_continue_additive == false)
		str_contact_faces = file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_contact_faces.txt";
	else
		str_contact_faces = file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock_contact_faces.txt";

	Visual vis_2;
	cv::Point3d input_ori(vector_after.x(), vector_after.y(), vector_after.z());
	if (judge_continue_additive == false)
		vis_2.generateArrows(input_ori, file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue));
	else
		vis_2.generateArrows(input_ori, file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(id_continue) + "_subblock");

	ofstream ofile_contact_faces(str_contact_faces);
	ofile_contact_faces << id_contact_faces.size() << endl;
	for (int i = 0; i < id_contact_faces.size(); i++)
		ofile_contact_faces << id_contact_faces[i] << endl;
	cout << "d" << endl;
	/////////////////////////
	return;
}

void HybridManufacturing::subtractive_accessibility_decomposition(
	vector<TRiangle> need_detect_triangle,
	int height_of_beam_search,
	int cont_number_of_queue,
	cutter cutting_tool,
	Slicer_2 current_slicer)
{
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 27;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 23;
	cutting_tool.carriage_height = 33;

	/////show all accessible points in every orientation/////
	Eigen::MatrixXd V_2;
	Eigen::MatrixXi F_2;
	igl::readOBJ("ball.obj", V_2, F_2);
	const Slicer_2 slicer = current_slicer; //加个const看看有没有涉及更改
	/*slicer.load((file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(cont_number_of_queue) + ".obj").c_str());*/

	/////////////////Get sample points//////////////////
	vector<cv::Point3d> all_sample_points_in_triangles(need_detect_triangle.size());
	for (int i = 0; i < need_detect_triangle.size(); i++) {

		cv::Point3d V_1(slicer.positions[need_detect_triangle[i][0]][0], slicer.positions[need_detect_triangle[i][0]][1], slicer.positions[need_detect_triangle[i][0]][2]);
		cv::Point3d V_2(slicer.positions[need_detect_triangle[i][1]][0], slicer.positions[need_detect_triangle[i][1]][1], slicer.positions[need_detect_triangle[i][1]][2]);
		cv::Point3d V_3(slicer.positions[need_detect_triangle[i][2]][0], slicer.positions[need_detect_triangle[i][2]][1], slicer.positions[need_detect_triangle[i][2]][2]);
		double a = distance3d(V_1, V_2);
		double b = distance3d(V_1, V_3);
		double c = distance3d(V_2, V_3);
		cv::Point3d V_incentre(((a * V_1.x + b * V_2.x + c * V_3.x) / (a + b + c)), ((a * V_1.y + b * V_2.y + c * V_3.y) / (a + b + c)), ((a * V_1.z + b * V_2.z + c * V_3.z) / (a + b + c)));
		all_sample_points_in_triangles[i] = V_incentre;
	}
	///////////////////////////////////////////////////
	cout << "% " << slicer.triangles.size() << endl;
	cout << "% " << need_detect_triangle.size() << endl;
	vector<vector<int>> accessible_ori_of_need_detect_V(need_detect_triangle.size(), vector<int>(sampling_subtractive.sample_points.size(), 10000000));
	vector<vector<Eigen::MatrixXd>> vis_points(1);
	vector<vector<vector<Eigen::Vector3d>>> vis_lines(1);


	cout << "%sampling_subtractive.sample_points.size() " << sampling_subtractive.sample_points.size() << endl;
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {
		//先去除水平面以下的方向  //注意只有最底下的block需要限制!!!!!!!!!!!!!!
		/*if (sampling_subtractive.sample_points[ori].z < 0.2)
			continue;*/
		vector<std::vector<Eigen::MatrixXd>> temp_new_V_remain(slicer.triangles.size());
		for (int i = 0; i < slicer.triangles.size(); i++) {
			temp_new_V_remain[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_new_V_remain[i][k].resize(3, 1);
				temp_new_V_remain[i][k](0, 0) = slicer.positions[slicer.triangles[i][k]][0];
				temp_new_V_remain[i][k](1, 0) = slicer.positions[slicer.triangles[i][k]][1];
				temp_new_V_remain[i][k](2, 0) = slicer.positions[slicer.triangles[i][k]][2];
			}
		}

		std::vector<Eigen::MatrixXd> temp_V_need_detect(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect[i].resize(3, 1);
			temp_V_need_detect[i](0, 0) = all_sample_points_in_triangles[i].x;
			temp_V_need_detect[i](1, 0) = all_sample_points_in_triangles[i].y;
			temp_V_need_detect[i](2, 0) = all_sample_points_in_triangles[i].z;
		}

		vector<std::vector<Eigen::MatrixXd>> temp_V_need_detect_triangle(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect_triangle[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_V_need_detect_triangle[i][k].resize(3, 1);
				temp_V_need_detect_triangle[i][k](0, 0) = slicer.positions[need_detect_triangle[i][k]][0];
				temp_V_need_detect_triangle[i][k](1, 0) = slicer.positions[need_detect_triangle[i][k]][1];
				temp_V_need_detect_triangle[i][k](2, 0) = slicer.positions[need_detect_triangle[i][k]][2];
			}
		}

		Eigen::Matrix3d rotMatrix;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < temp_new_V_remain.size(); i++)
			for (int j = 0; j < temp_new_V_remain[i].size(); j++)
				temp_new_V_remain[i][j] = rotMatrix.inverse() * temp_new_V_remain[i][j];
		for (int i = 0; i < temp_V_need_detect.size(); i++)
			temp_V_need_detect[i] = rotMatrix.inverse() * temp_V_need_detect[i];
		for (int i = 0; i < temp_V_need_detect_triangle.size(); i++)
			for (int j = 0; j < temp_V_need_detect_triangle[i].size(); j++)
				temp_V_need_detect_triangle[i][j] = rotMatrix.inverse() * temp_V_need_detect_triangle[i][j];

		vector<Eigen::Vector3d> all_normals_of_need_detect_triangle;
		all_normals_of_need_detect_triangle.clear();
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			cv::Point3d V_1(temp_V_need_detect_triangle[i][0](0, 0), temp_V_need_detect_triangle[i][0](1, 0), temp_V_need_detect_triangle[i][0](2, 0));
			cv::Point3d V_2(temp_V_need_detect_triangle[i][1](0, 0), temp_V_need_detect_triangle[i][1](1, 0), temp_V_need_detect_triangle[i][1](2, 0));
			cv::Point3d V_3(temp_V_need_detect_triangle[i][2](0, 0), temp_V_need_detect_triangle[i][2](1, 0), temp_V_need_detect_triangle[i][2](2, 0));
			double na = (V_2.y - V_1.y) * (V_3.z - V_1.z) - (V_2.z - V_1.z) * (V_3.y - V_1.y);
			double nb = (V_2.z - V_1.z) * (V_3.x - V_1.x) - (V_2.x - V_1.x) * (V_3.z - V_1.z);
			double nc = (V_2.x - V_1.x) * (V_3.y - V_1.y) - (V_2.y - V_1.y) * (V_3.x - V_1.x);
			Eigen::Vector3d vn(na, nb, nc);
			vn.normalize();
			all_normals_of_need_detect_triangle.push_back(vn);
			//cout << vn.x() << " "<<vn.y()<<" "<<vn.z() << endl;
		}

		vector<double> max_z_of_triangles(temp_new_V_remain.size());
		for (int i = 0; i < temp_new_V_remain.size(); i++) {
			max_z_of_triangles[i] = MIN_D;
			for (int j = 0; j < 3; j++)
				max_z_of_triangles[i] = max(max_z_of_triangles[i], temp_new_V_remain[i][j](2, 0));
		}

		///////////////////collision detection////////////////////////
		PrepareToolForCollision(cutting_tool);

		for (int i = 0; i < temp_V_need_detect.size(); i++) {
			bool jud_collision = false;
			Eigen::Vector3d center_point;

			center_point.x() = temp_V_need_detect[i](0, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].x();
			center_point.y() = temp_V_need_detect[i](1, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].y();
			center_point.z() = temp_V_need_detect[i](2, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].z();

			for (int ii = 0; ii < temp_new_V_remain.size(); ii++) {

				if (CheckToolCollisionWithCell(center_point, temp_new_V_remain[ii], max_z_of_triangles[ii], cutting_tool, 30.0, 3.0)) {
					jud_collision = true;
					break;
				}
			}

			if (jud_collision == false) {
				accessible_ori_of_need_detect_V[i][ori] = 0;
			}
		}


		//visualize
		/*rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorAfter,vectorBefore).toRotationMatrix();
		for (int i = 0; i < vis_points.size(); i++)
			vis_points[i] = rotMatrix.inverse() * vis_points[i];
		ofstream all_balls(".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + ".obj");
		for (int i = 0; i < vis_points.size(); i++) {
			for (int j = 0; j < V_2.rows(); j++)
				all_balls << "v " << V_2(j, 0) + vis_points[i](0, 0) << " " << V_2(j, 1) + vis_points[i](1, 0) << " " << V_2(j, 2) + vis_points[i](2, 0) << " 0.9" << " 0.05" << " 0.05" << endl;
			for (int j = 0; j < F_2.rows(); j++)
				all_balls << "f " << F_2(j, 0) + i * V_2.rows() + 1 << " " << F_2(j, 1) + i * V_2.rows() + 1 << " " << F_2(j, 2) + i * V_2.rows() + 1 << endl;
		}
		all_balls.close();
		Visual vis;
		vis.generateModelForRendering_7(sampling_subtractive.sample_points[ori], ".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + "_orientation.obj");*/
	}


	//////////////////////////////////graph cut////////////////////////////////////////
	//不应该出现cont_ori == 0的情况
	int cont_revise = 0;
	for (int i = 0; i < accessible_ori_of_need_detect_V.size(); i++) {
		int cont_ori = 0;
		for (int j = 0; j < accessible_ori_of_need_detect_V[i].size(); j++) {
			if (accessible_ori_of_need_detect_V[i][j] == 0)
				cont_ori++;
		}
		//暂时先强制任意方向可达
		if (cont_ori == 0) {
			cont_revise++;
			//cout << "*** " << need_detect_triangle[i][0] << endl;
			for (int j = 0; j < accessible_ori_of_need_detect_V[i].size(); j++)
				accessible_ori_of_need_detect_V[i][j] = 0;
		}
		/*if (accessible_ori_of_need_detect_V[i][123] != 0)
			cout << "no" << endl;*/
	}
	cout << cont_revise << endl;


	vector<std::vector<int>> pixels_relations;
	vector<int> length_edges;
	for (int i = 0; i < need_detect_triangle.size(); i++)
		for (int j = 0; j < need_detect_triangle.size(); j++) {
			bool jud_adjacent = false;
			vector<int> temp_edge(2);
			for (int ii = 0; ii < 3; ii++) {
				for (int jj = 0; jj < 3; jj++) {
					if (need_detect_triangle[i][ii] == need_detect_triangle[j][jj]) {
						/*temp_edge[0] = i;
						temp_edge[1] = j;
						pixels_relations.push_back(temp_edge);
						break;*/
						if (need_detect_triangle[i][(ii + 1) % 3] == need_detect_triangle[j][(jj + 1) % 3] || need_detect_triangle[i][(ii + 1) % 3] == need_detect_triangle[j][(jj + 2) % 3]) {
							temp_edge[0] = i;
							temp_edge[1] = j;
							pixels_relations.push_back(temp_edge);
							float temp_length = distanceVec3(slicer.positions[need_detect_triangle[i][ii]], slicer.positions[need_detect_triangle[i][(ii + 1) % 3]]) * 100;
							length_edges.push_back(temp_length);
							jud_adjacent = true;
							break;
						}
						else if (need_detect_triangle[i][(ii + 2) % 3] == need_detect_triangle[j][(jj + 1) % 3] || need_detect_triangle[i][(ii + 2) % 3] == need_detect_triangle[j][(jj + 2) % 3]) {
							temp_edge[0] = i;
							temp_edge[1] = j;
							pixels_relations.push_back(temp_edge);
							float temp_length = distanceVec3(slicer.positions[need_detect_triangle[i][ii]], slicer.positions[need_detect_triangle[i][(ii + 2) % 3]]) * 100;
							length_edges.push_back(temp_length);
							jud_adjacent = true;
							break;
						}
						else {
							jud_adjacent = true;
							break;
						}
					}
				}
				if (jud_adjacent == true)
					break;
			}
			/*if (jud_adjacent == false) {
				temp_edge[0] = i;
				temp_edge[1] = j;
				pixels_relations.push_back(temp_edge);
				length_edges.push_back(100000);
			}*/
		}


	cout << endl << "Graph Cut................" << endl;
	cout << "% " << need_detect_triangle.size() << endl;
	cout << "% " << sampling_subtractive.sample_points.size() << endl;
	vector<int> result = GeneralGraph_DArraySArraySpatVarying(need_detect_triangle.size(), sampling_subtractive.sample_points.size(), accessible_ori_of_need_detect_V, pixels_relations, length_edges);

	for (int i = 0; i < need_detect_triangle.size(); i++)
		for (int j = i + 1; j < need_detect_triangle.size(); j++)
			if (result[i] > result[j]) {
				swap(result[i], result[j]);
				swap(need_detect_triangle[i], need_detect_triangle[j]);
			}

	//cout << "normal:" << result[0]<<endl;
	ofstream ofile(".\\vis\\normal_of_pathches.txt");
	vector<vector<vasco::core::Tri3>> vis_triangles(1);
	vector<vasco::core::Vec3> vis_positions = slicer.positions;

	int cont_patch = 0;
	vector<Eigen::Vector3d> points_in_cell;
	vector<Eigen::Vector3d> normals;
	points_in_cell.push_back(Eigen::Vector3d(slicer.positions[need_detect_triangle[0][0]][0], slicer.positions[need_detect_triangle[0][0]][1], slicer.positions[need_detect_triangle[0][0]][2])); //point_in_cell初始放入第一个点
	normals.push_back(sampling_subtractive.sample_points[result[0]]);
	ofile << normals[0].x() << " " << normals[0].y() << " " << normals[0].z() << endl;


	for (int i = 0; i < need_detect_triangle.size(); i++) {
		if (i != 0 && result[i - 1] < result[i]) {
			cont_patch++;
			vector<Eigen::MatrixXd> temp_vecc(0);
			vis_points.push_back(temp_vecc);
			vector<vector<Eigen::Vector3d>> temp_vec(0);
			vis_lines.push_back(temp_vec);
			vector<vasco::core::Tri3> temp_tri(0);
			vis_triangles.push_back(temp_tri);
			points_in_cell.push_back(Eigen::Vector3d(slicer.positions[need_detect_triangle[i][0]][0], slicer.positions[need_detect_triangle[i][0]][1], slicer.positions[need_detect_triangle[i][0]][2]));
			normals.push_back(sampling_subtractive.sample_points[result[i]]);
			ofile << normals[normals.size() - 1].x() << " " << normals[normals.size() - 1].y() << " " << normals[normals.size() - 1].z() << endl;
		}

		//Eigen::MatrixXd temp_point(3, 1);
		//temp_point(0, 0) = V(index_V_need_detect[i], 0);
		//temp_point(1, 0) = V(index_V_need_detect[i], 1);
		//temp_point(2, 0) = V(index_V_need_detect[i], 2);
		////cout << vis_points[cont_patch].size() << endl;
		//vis_points[cont_patch].push_back(temp_point);
		//cout << "t" << endl;
		vector<Eigen::Vector3d> temp_vec;

		temp_vec.clear();
		vasco::core::Tri3 temp_tri;
		//cout << index_V_need_detect[i] << endl;
		for (int j = 0; j < 3; j++) {
			//cout << "f" << endl;
			Eigen::Vector3d temp_vec_2(slicer.positions[need_detect_triangle[i][j]][0], slicer.positions[need_detect_triangle[i][j]][1], slicer.positions[need_detect_triangle[i][j]][2]);
			//cout << "g" << endl;
			temp_vec.push_back(temp_vec_2);
			temp_tri[j] = need_detect_triangle[i][j];
		}
		//cout << "d" << endl;
		//cout << vis_lines.size() << " " << cont_patch << endl;
		vis_lines[cont_patch].push_back(temp_vec);
		vis_triangles[cont_patch].push_back(temp_tri);
		//cout << vis_lines[cont_patch].size() << endl;
		//cout << "e" << endl;
	}

	//visualize	
	/*ofstream all_balls(".\\vis\\coral_subtractive_decompose.obj");
	int cont_v = 0;
	for (int i = 0; i < vis_points.size(); i++) {
		double r = rand() / double(RAND_MAX);
		double g = rand() / double(RAND_MAX);
		double b = rand() / double(RAND_MAX);
		for (int j = 0; j < vis_points[i].size(); j++) {
			for (int k = 0; k < V_2.rows(); k++)
				all_balls << "v " << V_2(k, 0) + vis_points[i][j](0, 0) << " " << V_2(k, 1) + vis_points[i][j](1, 0) << " " << V_2(k, 2) + vis_points[i][j](2, 0) << " " << r << " " << g << " " << b << endl;
			for (int k = 0; k < F_2.rows(); k++)
				all_balls << "f " << F_2(k, 0) + cont_v * V_2.rows() + 1 << " " << F_2(k, 1) + cont_v * V_2.rows() + 1 << " " << F_2(k, 2) + cont_v * V_2.rows() + 1 << endl;
			cont_v++;
		}
	}
	all_balls.close();*/
	for (int t = 0; t < vis_lines.size(); t++) {
		std::ofstream dstream(".\\vis\\patch-" + to_string(height_of_beam_search) + "_" + to_string(t) + ".stl");
		if (!dstream.is_open()) {
			std::cout << "can not open " << std::endl;
			return;
		}
		dstream << "solid STL generated by MeshLab" << std::endl;
		for (int i = 0; i < vis_lines[t].size(); i++) {
			dstream << "  facet normal " << "0 0 0" << std::endl;
			dstream << "    outer loop" << std::endl;
			for (int j = 0; j < 3; j++) {
				dstream << "      vertex  " << vis_lines[t][i][j][0] << " " << vis_lines[t][i][j][1] << " " << vis_lines[t][i][j][2] << std::endl;
			}
			dstream << "    endloop" << std::endl;
			dstream << "  endfacet" << std::endl;
		}
		dstream << "endsolid vcg" << std::endl;
		dstream.close();
	}
	for (int t = 0; t < vis_triangles.size(); t++) {
		std::string mesh_name = "subtractive_patch_" + to_string(height_of_beam_search) + "_" + to_string(t);
		polyscope::registerSurfaceMesh(mesh_name, vis_positions, vis_triangles[t]);
	}
	polyscope::show();
	ofile.close();


	cout << cont_patch << endl;
	vector<vector<double>> colors(vis_lines.size(), vector<double>(3));
	Visual vis;
	vis.generateModelForRendering_10(vis_lines, colors, ".\\vis\\coral_subtractive_patch_voronoi_cell.obj");
	Visual vis_2;
	vis_2.generateModelForRendering_11(points_in_cell, normals, colors, ".\\vis\\coral_subtractive_patch_normal.obj");
	//Visual vis;
	//vis.generateModelForRendering_7(sampling_subtractive.sample_points[ori], ".\\vis\\coral_accessible_points_in_ori-" + to_string(ori) + "_orientation.obj");
	//////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////
}

void HybridManufacturing::subtractive_accessibility_decomposition_within_2_blocks(int height_of_beam_search, cutter cutting_tool)
{
	sampling_subtractive.OrientationSamplePoints();	//sampling_subtractive生成球面采样点
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 27;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 23;
	cutting_tool.carriage_height = 33;

	/////show all accessible points in every orientation/////
	Eigen::MatrixXd V_2;
	Eigen::MatrixXi F_2;
	igl::readOBJ("ball.obj", V_2, F_2);

	Slicer_2 slicer_load_current_patch;

	Slicer_2 slicer_load_next_patch;

	// 读取并拼接 .\vis\block_patch-1_.obj ~ .\vis\block_patch-height_of_beam_search_.obj
	std::vector<int> merged_vertex_source_patch_id;
	std::vector<int> merged_face_source_patch_id;

	Slicer_2 slicer_load_merged_patch = MergeBlockPatchesWithDedup(
		height_of_beam_search,
		merged_vertex_source_patch_id,
		merged_face_source_patch_id,
		1e-6);
	slicer_load_merged_patch.save(".\\vis\\block_patch_merged_removedup.obj");
	std::cout << "[Info] merged patch mesh: V=" << slicer_load_merged_patch.positions.size()
		<< ", F=" << slicer_load_merged_patch.triangles.size() << std::endl;

	// 计算 merged patch 上每个三角面在每个采样方向下的最大碰撞 patch_index
	const auto merged_face_min_collision_patch = EvaluateMergedPatchToolCollision(
		slicer_load_merged_patch,
		merged_face_source_patch_id,
		cutting_tool);

	std::cout << "[Info] merged-face min-collision-patch matrix computed: faces="
		<< merged_face_min_collision_patch.size() << ", orientations="
		<< sampling_subtractive.sample_points.size() << std::endl;

	// ---------------- graph cut on merged patch ----------------
	{
		const int face_count = static_cast<int>(slicer_load_merged_patch.triangles.size());
		const int ori_count = static_cast<int>(sampling_subtractive.sample_points.size());
		const int block_count = height_of_beam_search;
		const int num_labels = block_count * ori_count;
		const int INF_COST = 10000000;

		if (face_count == 0 || ori_count == 0 || block_count <= 0) {
			std::cout << "[Warn] skip graph cut: invalid sizes. face_count=" << face_count
				<< ", ori_count=" << ori_count << ", block_count=" << block_count << std::endl;
		}
		else {
			auto encode_label = [ori_count](int patch_id, int ori_id) -> int {
				// patch_id: 1..block_count, ori_id: 0..ori_count-1
				return (patch_id - 1) * ori_count + ori_id;
				};

			auto decode_label = [ori_count](int label_id, int& patch_id, int& ori_id) {
				patch_id = label_id / ori_count + 1;
				ori_id = label_id % ori_count;
				};

			// data cost: face x label
			std::vector<std::vector<int>> data_value(face_count, std::vector<int>(num_labels, INF_COST));

			int warn_a_lt_b = 0;
			int warn_no_feasible = 0;

			for (int i = 0; i < face_count; ++i) {
				int b = 1;
				if (i < static_cast<int>(merged_face_source_patch_id.size())) {
					b = merged_face_source_patch_id[i];
				}
				else {
					std::cout << "[what?]i >= static_cast<int>(merged_face_source_patch_id.size())" << std::endl;
				}
				b = std::max(1, std::min(block_count, b));

				int feasible_cnt = 0;

				for (int ori = 0; ori < ori_count; ++ori) {
					int a = merged_face_min_collision_patch[i][ori];
					if (a == -1) {
						a = 1;
					}

					a = std::max(1, std::min(a, block_count));

					// 按 [a..b] 添加标签，若 a>b 不添加并提示
					if (a > b) {
						++warn_a_lt_b;
						continue;
					}

					for (int p = a; p <= b; ++p) {
						const int lid = encode_label(p, ori);
						data_value[i][lid] = 0;
						++feasible_cnt;
					}
				}

				// 防御：若该面没有任何可行label，给一个保底label避免graph cut退化
				if (feasible_cnt == 0) {
					++warn_no_feasible;
					data_value[i][encode_label(b, 30)] = 0;
				}
			}

			if (warn_a_lt_b > 0) {
				std::cout << "[Warn] graph-cut label range invalid (a<b) count = " << warn_a_lt_b << std::endl;
			}
			if (warn_no_feasible > 0) {
				std::cout << "[Warn] faces with no feasible label = " << warn_no_feasible
					<< " (fallback label assigned)." << std::endl;
			}

			// 建立面邻接图（共享边）
			std::vector<std::vector<int>> pixels_relations;
			std::vector<int> length_edges;

			std::map<std::pair<int, int>, int> edge_owner; // edge -> face id
			edge_owner.clear();

			for (int i = 0; i < face_count; ++i) {
				const auto& tri = slicer_load_merged_patch.triangles[i];
				for (int e = 0; e < 3; ++e) {
					int u = tri[e];
					int v = tri[(e + 1) % 3];
					if (u > v) std::swap(u, v);

					const std::pair<int, int> key(u, v);
					auto it = edge_owner.find(key);
					if (it == edge_owner.end()) {
						edge_owner.insert({ key, i });
					}
					else {
						const int j = it->second;
						if (j != i) {
							pixels_relations.push_back({ j, i });

							const auto& p1 = slicer_load_merged_patch.positions[u];
							const auto& p2 = slicer_load_merged_patch.positions[v];
							const int w = std::max(1, static_cast<int>(distanceVec3(p1, p2) * 100.0));
							length_edges.push_back(w);
						}
					}
				}
			}

			std::cout << "[Info] graph-cut nodes=" << face_count
				<< ", labels=" << num_labels
				<< ", edges=" << pixels_relations.size() << std::endl;

			// graph cut
			std::vector<int> result_labels = GeneralGraph_DArraySArraySpatVarying(
				face_count,
				num_labels,
				data_value,
				pixels_relations,
				length_edges);

			// 可视化：按最终label所属patch着色（离散固定色）
			const std::string gc_obj_file = ".\\vis\\merged_patch_graphcut_label.obj";
			std::ofstream ofs(gc_obj_file);
			if (!ofs.is_open()) {
				std::cout << "[Warn] cannot open file for writing: " << gc_obj_file << std::endl;
			}
			else {
				auto color_from_patch = [](int ori_id) -> std::array<double, 3> {
					static const std::array<std::array<double, 3>, 32> palette = { {
						{0.894, 0.102, 0.110}, {0.216, 0.494, 0.722}, {0.302, 0.686, 0.290}, {0.596, 0.306, 0.639},
						{1.000, 0.498, 0.000}, {1.000, 1.000, 0.200}, {0.651, 0.337, 0.157}, {0.969, 0.506, 0.749},
						{0.600, 0.600, 0.600}, {0.121, 0.466, 0.705}, {0.682, 0.780, 0.909}, {1.000, 0.733, 0.470},
						{0.172, 0.627, 0.172}, {0.839, 0.153, 0.157}, {0.580, 0.404, 0.741}, {0.549, 0.337, 0.294},
						{0.890, 0.467, 0.761}, {0.498, 0.498, 0.498}, {0.737, 0.741, 0.133}, {0.090, 0.745, 0.811},
						{0.400, 0.760, 0.647}, {0.988, 0.553, 0.384}, {0.553, 0.627, 0.796}, {0.906, 0.541, 0.765},
						{0.651, 0.847, 0.329}, {1.000, 0.851, 0.184}, {0.898, 0.768, 0.580}, {0.702, 0.702, 0.702},
						{0.984, 0.603, 0.600}, {0.800, 0.922, 0.773}, {0.871, 0.796, 0.894}, {0.996, 0.851, 0.651}
					} };

					int idx = ori_id;
					if (idx < static_cast<int>(palette.size())) {
						return palette[static_cast<size_t>(idx)];
					}

					return palette[idx % palette.size()];

					};

				// 同时输出 face -> (patch,ori)
				std::ofstream lfs(".\\vis\\merged_patch_graphcut_label.txt");
				lfs << "face_id patch_id ori_id label_id\n";

				int v_count = 0;
				int skipped_faces = 0;
				int invalid_vertices = 0;

				// 按 patch_id 收集面，供 polyscope 分 mesh 可视化
				std::map<int, std::vector<vasco::core::Tri3>> patch_to_faces;

				for (int i = 0; i < face_count; ++i) {
					int patch_id = 1, ori_id = 0;
					decode_label(result_labels[i], patch_id, ori_id);

					// 修正：按 patch_id 取颜色
					const auto color = color_from_patch(result_labels[i]);
					const auto& tri = slicer_load_merged_patch.triangles[i];

					// tri 索引与坐标合法性检查
					bool face_ok = true;
					for (int k = 0; k < 3; ++k) {
						const int vid = tri[k];
						if (vid < 0 || vid >= static_cast<int>(slicer_load_merged_patch.positions.size())) {
							face_ok = false;
							break;
						}
						const auto& p = slicer_load_merged_patch.positions[vid];
						if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
							face_ok = false;
							++invalid_vertices;
							break;
						}
					}

					if (!face_ok) {
						++skipped_faces;
						continue;
					}

					for (int k = 0; k < 3; ++k) {
						const auto& p = slicer_load_merged_patch.positions[tri[k]];
						ofs << "v " << p[0] << " " << p[1] << " " << p[2] << " "
							<< color[0] << " " << color[1] << " " << color[2] << "\n";
					}
					ofs << "f " << (v_count + 1) << " " << (v_count + 2) << " " << (v_count + 3) << "\n";
					v_count += 3;

					lfs << i << " " << patch_id << " " << ori_id << " " << result_labels[i] << "\n";

					// 保存到 patch 子网格
					patch_to_faces[result_labels[i]].push_back(tri);
				}

				if (skipped_faces > 0) {
					std::cout << "[Warn] skipped faces while writing OBJ: " << skipped_faces
						<< ", invalid_vertices=" << invalid_vertices << std::endl;
				}

				// polyscope: 按 patch_id 拆分 mesh，并设置固定颜色
				std::vector<Eigen::Vector3d> patch_arrow_points;
				std::vector<Eigen::Vector3d> patch_arrow_dirs;
				std::vector<std::vector<double>> patch_arrow_colors;

				Eigen::Vector3d bb_min(MAX_D, MAX_D, MAX_D), bb_max(MIN_D, MIN_D, MIN_D);
				for (const auto& p : slicer_load_merged_patch.positions) {
					bb_min.x() = std::min(bb_min.x(), p[0]);
					bb_min.y() = std::min(bb_min.y(), p[1]);
					bb_min.z() = std::min(bb_min.z(), p[2]);
					bb_max.x() = std::max(bb_max.x(), p[0]);
					bb_max.y() = std::max(bb_max.y(), p[1]);
					bb_max.z() = std::max(bb_max.z(), p[2]);
				}
				const double model_diag = (bb_max - bb_min).norm();
				const double arrow_start_offset = std::max(1.0, 0.06 * model_diag);
				const double arrow_length = std::max(2.0, 0.12 * model_diag);

				for (const auto& kv : patch_to_faces) {
					const int patch_id = kv.first;
					const auto& patch_faces = kv.second;
					if (patch_faces.empty()) continue;
					int de_patch;
					int ori_id = 0;
					decode_label(patch_id, de_patch, ori_id);
					const std::string mesh_name =
						"merged_patch_graphcut_patch_" + std::to_string(height_of_beam_search) + "_" + std::to_string(de_patch) + "_" + std::to_string(ori_id);

					auto* ps_mesh = polyscope::registerSurfaceMesh(
						mesh_name,
						slicer_load_merged_patch.positions,
						patch_faces);

					const auto c = color_from_patch(ori_id);
					ps_mesh->setSurfaceColor({ (float)c[0], (float)c[1], (float)c[2] });

					std::cout << "[debug]mesh_name = " << mesh_name << ", patch_id=" << patch_id << ", de_patch=" << de_patch << ", ori_id=" << ori_id
						<< ", c= " << c[0] << "," << c[1] << "," << c[2] << std::endl;

					Eigen::Vector3d patch_center(0.0, 0.0, 0.0);
					int center_count = 0;
					for (const auto& tri : patch_faces) {
						for (int k = 0; k < 3; ++k) {
							const int vid = tri[k];
							if (vid >= 0 && vid < static_cast<int>(slicer_load_merged_patch.positions.size())) {
								const auto& p = slicer_load_merged_patch.positions[vid];
								patch_center.x() += p[0];
								patch_center.y() += p[1];
								patch_center.z() += p[2];
								++center_count;
							}
						}
					}
					if (center_count > 0) {
						patch_center /= static_cast<double>(center_count);

						Eigen::Vector3d ori_dir(0.0, 0.0, 1.0);
						if (ori_id >= 0 && ori_id < static_cast<int>(sampling_subtractive.sample_points.size())) {
							ori_dir = sampling_subtractive.sample_points[ori_id];
						}
						Eigen::Vector3d ori_vec(ori_dir);
						double n = ori_vec.norm();
						if (n < 1e-12) {
							ori_vec = Eigen::Vector3d(0.0, 0.0, 1.0);
						}
						else {
							ori_vec /= n;
						}

						const Eigen::Vector3d arrow_start = patch_center + ori_vec * arrow_start_offset;
						const Eigen::Vector3d arrow_end = arrow_start + ori_vec * arrow_length;

						patch_arrow_points.push_back(arrow_start);
						patch_arrow_dirs.push_back(ori_vec);

						patch_arrow_colors.push_back({ c[0], c[1], c[2] });

						const std::string arrow_name =
							"merged_patch_graphcut_arrow_" + std::to_string(height_of_beam_search) + "_" + std::to_string(de_patch) + "_" + std::to_string(ori_id);

						Eigen::Vector3d ref_axis = (std::abs(ori_vec.z()) < 0.9)
							? Eigen::Vector3d(0.0, 0.0, 1.0)
							: Eigen::Vector3d(0.0, 1.0, 0.0);
						Eigen::Vector3d perp1 = ori_vec.cross(ref_axis);
						double perp1_norm = perp1.norm();
						if (perp1_norm < 1e-12) {
							perp1 = Eigen::Vector3d(1.0, 0.0, 0.0);
						}
						else {
							perp1 /= perp1_norm;
						}

						const double shaft_w = std::max(0.2, arrow_length * 0.08);
						const double head_len = arrow_length * 0.35;
						const double head_w = shaft_w * 2.0;

						Eigen::Vector3d m = arrow_end - ori_vec * head_len;
						Eigen::Vector3d s0 = arrow_start + perp1 * shaft_w;
						Eigen::Vector3d s1 = arrow_start - perp1 * shaft_w;
						Eigen::Vector3d m0 = m + perp1 * shaft_w;
						Eigen::Vector3d m1 = m - perp1 * shaft_w;
						Eigen::Vector3d h0 = m + perp1 * head_w;
						Eigen::Vector3d h1 = m - perp1 * head_w;

						std::vector<vasco::core::Vec3> arrow_pos = {
							{ s0.x(), s0.y(), s0.z() },
							{ s1.x(), s1.y(), s1.z() },
							{ m1.x(), m1.y(), m1.z() },
							{ m0.x(), m0.y(), m0.z() },
							{ h0.x(), h0.y(), h0.z() },
							{ h1.x(), h1.y(), h1.z() },
							{ arrow_end.x(), arrow_end.y(), arrow_end.z() }
						};
						std::vector<vasco::core::Tri3> arrow_tri = {
							{ 0, 1, 2 }, { 0, 2, 3 },
							{ 4, 5, 6 }, { 3, 4, 6 }, { 2, 6, 5 }
						};

						auto* ps_arrow = polyscope::registerSurfaceMesh(arrow_name, arrow_pos, arrow_tri);
						ps_arrow->setSurfaceColor({ (float)c[0], (float)c[1], (float)c[2] });
					}
				}

				if (!patch_arrow_points.empty()) {
					Visual vis_patch_ori;
					vis_patch_ori.generateModelForRendering_11(
						patch_arrow_points,
						patch_arrow_dirs,
						patch_arrow_colors,
						".\\vis\\merged_patch_graphcut_patch_ori_arrows.obj");
				}
				polyscope::show();

				std::cout << "[Info] wrote graph-cut colored OBJ: " << gc_obj_file << std::endl;
				std::cout << "[Info] wrote graph-cut labels txt: .\\vis\\merged_patch_graphcut_label.txt" << std::endl;
			}
		}
	}

}

vector<vector<int>> HybridManufacturing::getAccessOri(const Slicer_2& slicer, Slicer_2& slicer_load_patch, vector<vasco::core::Vec3>& all_sample_points_in_triangles, cutter cutting_tool)
{
	auto& need_detect_triangle = slicer_load_patch.triangles;
	vector<vector<int>> accessible_ori_of_need_detect_V(need_detect_triangle.size(), vector<int>(sampling_subtractive.sample_points.size(), 10000000));
	cout << "%sampling_subtractive.sample_points.size() " << sampling_subtractive.sample_points.size() << endl;
	for (int ori = 0; ori < sampling_subtractive.sample_points.size(); ori++) {
		vector<std::vector<Eigen::MatrixXd>> temp_new_V_remain(slicer.triangles.size());
		for (int i = 0; i < slicer.triangles.size(); i++) {
			temp_new_V_remain[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_new_V_remain[i][k].resize(3, 1);
				temp_new_V_remain[i][k](0, 0) = slicer.positions[slicer.triangles[i][k]][0];
				temp_new_V_remain[i][k](1, 0) = slicer.positions[slicer.triangles[i][k]][1];
				temp_new_V_remain[i][k](2, 0) = slicer.positions[slicer.triangles[i][k]][2];
			}
		}

		std::vector<Eigen::MatrixXd> temp_V_need_detect(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect[i].resize(3, 1);
			temp_V_need_detect[i](0, 0) = all_sample_points_in_triangles[i][0];
			temp_V_need_detect[i](1, 0) = all_sample_points_in_triangles[i][1];
			temp_V_need_detect[i](2, 0) = all_sample_points_in_triangles[i][2];
		}

		vector<std::vector<Eigen::MatrixXd>> temp_V_need_detect_triangle(need_detect_triangle.size());
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			temp_V_need_detect_triangle[i].resize(3);
			for (int k = 0; k < 3; k++) {
				temp_V_need_detect_triangle[i][k].resize(3, 1);
				temp_V_need_detect_triangle[i][k](0, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][0];
				temp_V_need_detect_triangle[i][k](1, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][1];
				temp_V_need_detect_triangle[i][k](2, 0) = slicer_load_patch.positions[need_detect_triangle[i][k]][2];
			}
		}

		Eigen::Matrix3d rotMatrix;
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling_subtractive.sample_points[ori]);
		vectorAfter.normalize();
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		for (int i = 0; i < temp_new_V_remain.size(); i++)
			for (int j = 0; j < temp_new_V_remain[i].size(); j++)
				temp_new_V_remain[i][j] = rotMatrix.inverse() * temp_new_V_remain[i][j];
		for (int i = 0; i < temp_V_need_detect.size(); i++)
			temp_V_need_detect[i] = rotMatrix.inverse() * temp_V_need_detect[i];
		for (int i = 0; i < temp_V_need_detect_triangle.size(); i++)
			for (int j = 0; j < temp_V_need_detect_triangle[i].size(); j++)
				temp_V_need_detect_triangle[i][j] = rotMatrix.inverse() * temp_V_need_detect_triangle[i][j];

		vector<Eigen::Vector3d> all_normals_of_need_detect_triangle;
		all_normals_of_need_detect_triangle.clear();
		for (int i = 0; i < need_detect_triangle.size(); i++) {
			cv::Point3d V_1(temp_V_need_detect_triangle[i][0](0, 0), temp_V_need_detect_triangle[i][0](1, 0), temp_V_need_detect_triangle[i][0](2, 0));
			cv::Point3d V_2(temp_V_need_detect_triangle[i][1](0, 0), temp_V_need_detect_triangle[i][1](1, 0), temp_V_need_detect_triangle[i][1](2, 0));
			cv::Point3d V_3(temp_V_need_detect_triangle[i][2](0, 0), temp_V_need_detect_triangle[i][2](1, 0), temp_V_need_detect_triangle[i][2](2, 0));
			double na = (V_2.y - V_1.y) * (V_3.z - V_1.z) - (V_2.z - V_1.z) * (V_3.y - V_1.y);
			double nb = (V_2.z - V_1.z) * (V_3.x - V_1.x) - (V_2.x - V_1.x) * (V_3.z - V_1.z);
			double nc = (V_2.x - V_1.x) * (V_3.y - V_1.y) - (V_2.y - V_1.y) * (V_3.x - V_1.x);
			Eigen::Vector3d vn(na, nb, nc);
			vn.normalize();
			all_normals_of_need_detect_triangle.push_back(vn);
			//cout << vn.x() << " "<<vn.y()<<" "<<vn.z() << endl;
		}

		vector<double> max_z_of_triangles(temp_new_V_remain.size());
		for (int i = 0; i < temp_new_V_remain.size(); i++) {
			max_z_of_triangles[i] = MIN_D;
			for (int j = 0; j < 3; j++)
				max_z_of_triangles[i] = max(max_z_of_triangles[i], temp_new_V_remain[i][j](2, 0));
		}

		///////////////////collision detection////////////////////////
		PrepareToolForCollision(cutting_tool);

		for (int i = 0; i < temp_V_need_detect.size(); i++) {
			Eigen::Vector3d center_point;

			center_point.x() = temp_V_need_detect[i](0, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].x();
			center_point.y() = temp_V_need_detect[i](1, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].y();
			center_point.z() = temp_V_need_detect[i](2, 0) + (cutting_tool.cylinder_r) * all_normals_of_need_detect_triangle[i].z();

			bool jud_collision = false;

			for (int ii = 0; ii < temp_new_V_remain.size(); ii++) {
				if (CheckToolCollisionWithCell(center_point, temp_new_V_remain[ii], max_z_of_triangles[ii], cutting_tool, 30.0, 3.0)) {
					jud_collision = true;
					break;
				}
			}

			if (jud_collision == false) {
				//vis_points.push_back(temp_V_need_detect[i]);
				accessible_ori_of_need_detect_V[i][ori] = 0;
			}
			//else
				//cout << "% " << endl;
		}
	}
	return accessible_ori_of_need_detect_V;
}

void HybridManufacturing::outer_beam_search(nozzle the_nozzle, cutter cutting_tool)
{
	int W1 = 1;  //4
	clock_t start_time_total, end_time_total;
	clock_t start_time, end_time;
	clock_t start_time_2, end_time_2;
	clock_t start_time_3, end_time_3;
	clock_t start_time_4, end_time_4;
	clock_t start_time_5, end_time_5;
	clock_t start_time_6, end_time_6;
	clock_t start_time_7, end_time_7;
	clock_t start_time_8, end_time_8;
	float total_time_1 = 0, total_time_2 = 0, total_time_3 = 0, total_time_4 = 0;
	double sum_time = 0;
	double sum_time_2 = 0;
	double sum_time_3 = 0;
	double sum_time_4 = 0;
	double sum_time_5 = 0;
	double sum_time_6 = 0;
	int Sum_candidate_blocks = 0;
	int Sum_connected_components = 0;
	vector<double> evaluation_value(6, 0);
	cont_extra_additive_orientation = 0;

	start_time_total = clock();
	bool jud_outer_beam_search_terminate = false;
	vector<vector<vector<cv::Point3d>>> Tree_nodes;
	vector<vector<vector<cv::Point3d>>> Tree_nodes_contain;
	vector<vector<int>> Tree_nodes_cut_layers;
	vector<vector<int>> Tree_nodes_num_of_cut_layers_dependency_layer;
	vector<vector<Eigen::MatrixXd>> Tree_nodes_fragile_V;
	vector<double> Tree_nodes_larger_base;
	vector<vector<vector<area_S>>> Tree_nodes_ori_all_the_area_S;
	vector<vector<bool>> Tree_nodes_judge_S_be_searched;
	vector<bool> Tree_nodes_judge_continue;
	vector<int> Tree_nodes_continue_id;
	vector<bool> Tree_nodes_error;
	int* pre_tree_nodes = new int[500000];
	int* pre_index_of_nodes = new int[500000];
	vector<int> candidate_nodes;	//当前层生成的候选节点
	queue<int> last_step_nodes;
	vector<vector<Eigen::Vector3d>> save_ori;
	memset(pre_tree_nodes, -1, sizeof(pre_tree_nodes));
	vector<vector<cv::Point3d>> root_node;
	vector<int> root_node_2;
	double root_node_3;
	vector<Eigen::MatrixXd> root_node_4;
	Eigen::Vector3d root_ori;
	vector<vector<area_S>> root_node_5;
	vector<bool> root_node_6;
	vector<Eigen::Vector3d> root_node_7;
	root_node_7.push_back(root_ori);	//推入根节点的初始方向
	Tree_nodes.push_back(root_node);
	Tree_nodes_contain.push_back(root_node);
	Tree_nodes_cut_layers.push_back(root_node_2);
	Tree_nodes_num_of_cut_layers_dependency_layer.push_back(root_node_2);
	Tree_nodes_fragile_V.push_back(root_node_4);
	Tree_nodes_larger_base.push_back(root_node_3);
	Tree_nodes_ori_all_the_area_S.push_back(root_node_5);
	Tree_nodes_judge_S_be_searched.push_back(root_node_6);
	last_step_nodes.push(0);	//推入根节点，作为第一层搜索的初始节点
	Tree_nodes_judge_continue.push_back(false);
	Tree_nodes_continue_id.push_back(-1);
	Tree_nodes_error.push_back(false);
	save_ori.push_back(root_node_7);
	vector<vector<area_S>> ori_all_the_area_S = all_the_area_S;	//不可达点的all_the_area_S的拷贝
	ori_all_the_covering_points = all_the_covering_points;	//不可达点的all_the_covering_points的拷贝
	SAMPLE_ON_BALL sampling;
	int index_node = 0;
	sampling.OrientationSamplePoints_2();	//生成球面采样方向，作为增材打印方向候选
	//Visual Vis_ori;
	//Vis_ori.generateModelForRendering_6(sampling.sample_points, file_name);


	GetALLFragileVertex(sampling);

	Eigen::Matrix3d rotMatrix;
	vector<vector<int>> final_pathes_include_S;
	vector<vector<int>> final_pathes_include_sample_points;
	vector<vector<int>> final_pathes_include_covering_points;
	final_pathes_include_S.push_back(root_node_2);
	final_pathes_include_sample_points.push_back(root_node_2);
	final_pathes_include_covering_points.push_back(root_node_2);
	all_saved_mesh.resize(W1);

	int cont_number_of_queue = 0;	//记录当前处理到第几个节点，用于设定读取的文件名等
	int height_of_beam_search = 0;

	ofstream ofile_cont("E:\\Hybrid manufacturing\\HybridManufacturing\\HybridManufacturing\\models\\coral\\cont_red_number.txt");
	ofstream ofile_cont_size("E:\\Hybrid manufacturing\\HybridManufacturing\\HybridManufacturing\\models\\coral\\cont_size.txt");
	/*bool* flag_sample_point_used = new bool[sampling.sample_points.size()];
	for (int i = 0; i < sampling.sample_points.size(); i++)
		flag_sample_point_used[i] = false;*/
	cout << "--------------------- height of outer-beam search: 1 ---------------------";
	while (last_step_nodes.size() != 0 || candidate_nodes.size() != 0) {
		/*if (height_of_beam_search ==2)
			W1 = 5;
		else
			W1 = 1;*/
		start_time_2 = clock();
		vector<vector<bool>> is_fragile_V_2 = is_fragile_V;	//fragile信息的快照
		int now_last_node = last_step_nodes.front();
		/*if (Tree_nodes_error[now_last_node] == true) {
			last_step_nodes.pop();
			continue;
		}*/

		cout << endl << "Slicing and inner-beam-search......" << endl;
		//load blocks//
		const char* config_path = "config.ini";
		Katana::Instance().config.loadConfig(config_path);
		if (Tree_nodes_judge_continue[now_last_node] == false) {
			Katana::Instance().stl.loadStl((file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + ".stl").c_str());	//加载当前节点对应的模型
			cout << "T" << endl;
		}
		else {
			Katana::Instance().stl.loadStl((file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(pre_index_of_nodes[now_last_node]) + ".stl").c_str());
			cout << "U" << now_last_node << endl;
		}
		Katana::Instance().temp_vertices.resize(Katana::Instance().vertices.size());	//更新当前模型的顶点信息
		Katana::Instance().temp_vertices = Katana::Instance().vertices;
		end_time_2 = clock();
		sum_time_5 += double(end_time_2 - start_time_2) / CLOCKS_PER_SEC;
		//////////////

		start_time_5 = clock();
		//update subtractive dependency graph//
		vector<int> index_V_in_the_remaining_blocks;
		vector<int> S_in_block;
		vector<int> sample_points_in_block;	//记录不存在于当前块中的顶点索引
		vector<int> covering_points_in_block;	//记录被碰撞单元索引
		vector<bool> judge_S_be_searched; vector<bool> judge_covering_points_be_searched;
		judge_S_be_searched.clear(); judge_covering_points_be_searched.clear();
		judge_S_be_searched.resize(all_the_area_S.size());
		judge_covering_points_be_searched.resize(all_the_covering_points.size());
		for (int i = 0; i < all_the_area_S.size(); i++)
			judge_S_be_searched[i] = false;
		for (int i = 0; i < all_the_covering_points.size(); i++)
			judge_covering_points_be_searched[i] = false;
		S_in_block.clear(); sample_points_in_block.clear(); covering_points_in_block.clear(); index_V_in_the_remaining_blocks.clear();
		for (int i = 0; i < V.rows(); i++) {	//枚举原始输入模型网格V的所有顶点，判断该顶点是否仍存在于当前块中
			bool jud_still_exist = false;
			/*if (flag_sample_point_used[i])
				continue;*/
			for (int j = 0; j < Katana::Instance().vertices.size(); j++) {
				if (abs(V(i, 0) - Katana::Instance().vertices[j].x) <= 0.001 && abs(V(i, 1) - Katana::Instance().vertices[j].y) <= 0.001 && abs(V(i, 2) - Katana::Instance().vertices[j].z) <= 0.001) {
					index_V_in_the_remaining_blocks.push_back(i);	//如果判断该顶点仍存在于当前块中，则将其索引i加入index_V_in_the_remaining_blocks
					jud_still_exist = true;	//标记该顶点仍存在
					//flag_sample_point_used[i] = true;
					break;
				}
			}

			if (jud_still_exist == false) {	//如果该顶点不存在于当前块中
				for (int k = 0; k < is_fragile_V_2.size(); k++)
					is_fragile_V_2[k][i] = false;	//is_fragile_V_2中对应该顶点的信息全部置为false
				sample_points_in_block.push_back(i);
				for (int k = 0; k < map_S_and_vertex.size(); k++) {
					if (map_S_and_vertex[k] == i)
						S_in_block.push_back(k);	//将该顶点对应的不可达索引加入S_in_block
				}
				for (int k = 0; k < map_covering_points_and_vertex.size(); k++) {
					if (map_covering_points_and_vertex[k] == i)
						covering_points_in_block.push_back(k);	//将该顶点对应的被碰撞单元索引加入covering_points_in_block
				}
			}
		}

		/*while (pre_tree_nodes[now_last_node] != -1) {
			for (int i = 0; i < final_pathes_include_S[now_last_node].size(); i++) {
				S_in_block.push_back(final_pathes_include_S[now_last_node][i]);
			}
			for (int i = 0; i < final_pathes_include_sample_points[now_last_node].size(); i++) {
				sample_points_in_block.push_back(final_pathes_include_sample_points[now_last_node][i]);
			}
			now_last_node = pre_tree_nodes[now_last_node];
		}
		now_last_node = last_step_nodes.front();*/

		for (int i = 0; i < S_in_block.size(); i++)
			judge_S_be_searched[S_in_block[i]] = true;

		ori_num_points_of_ori_in_all_the_area_S.clear();
		ori_num_points_of_ori_in_all_the_area_S.resize(ori_all_the_area_S.size());	//为每个不可达点在ori_num_points_of_ori_in_all_the_area_S分配一个向量
		for (int i = 0; i < ori_all_the_area_S.size(); i++) {
			ori_num_points_of_ori_in_all_the_area_S[i].resize(sampling_subtractive.sample_points.size());	//每个不可达点记录在各个采样方向上对应的点数
			for (int j = 0; j < ori_num_points_of_ori_in_all_the_area_S[i].size(); j++)
				ori_num_points_of_ori_in_all_the_area_S[i][j] = 0;
		}
		for (int i = 0; i < ori_all_the_area_S.size(); i++) {	//枚举每个不可达点
			for (int itr = 0; itr < ori_all_the_area_S[i].size(); itr++) {	//枚举该不可达点对应的所有area_S
				ori_num_points_of_ori_in_all_the_area_S[i][ori_all_the_area_S[i][itr].oriId]++;	//统计该不可达点在各个采样方向上对应的点数
			}
		}
		for (int i = 0; i < sample_points_in_block.size(); i++) {  //枚举每个不存在于当前块中的顶点
			for (int j = 0; j < all_the_covering_points[map_covering_points_and_vertex_inv[sample_points_in_block[i]]].size(); j++) {
				int index = all_the_covering_points[map_covering_points_and_vertex_inv[sample_points_in_block[i]]][j].pointId;
				int ori = all_the_covering_points[map_covering_points_and_vertex_inv[sample_points_in_block[i]]][j].oriId;
				ori_num_points_of_ori_in_all_the_area_S[index][ori]--;	//更新该不可达点在该采样方向上对应的点数
			}
		}
		for (int i = 0; i < covering_points_in_block.size(); i++)
			judge_covering_points_be_searched[covering_points_in_block[i]] = true;	//标记需要更新的被碰撞单元

		if (Tree_nodes_judge_continue[now_last_node] == true) {
			Katana::Instance().stl.loadStl((file_name + "-" + to_string(height_of_beam_search) + "_" + to_string(cont_number_of_queue) + "_" + to_string(Tree_nodes_continue_id[now_last_node]) + "_subblock.stl").c_str());
			Katana::Instance().temp_vertices.resize(Katana::Instance().vertices.size());
			Katana::Instance().temp_vertices = Katana::Instance().vertices;
		}
		end_time_5 = clock();
		sum_time_3 += double(end_time_5 - start_time_5) / CLOCKS_PER_SEC;
		/////////////////////////////////////

		/*int end;
		if (height_of_beam_search != 4)
			end = sampling.sample_points.size();
		else
			end = 1;*/
		std::cout << sampling.sample_points.size() << endl;
		for (int ori = 0; ori < sampling.sample_points.size(); ori++) {	//枚举所有采样方向，作为增材分层方向
			printf("\r[%d%%]>", ori * 100 / (sampling.sample_points.size() - 1));
			for (int j = 1; j <= ori * 20 / sampling.sample_points.size(); j++)
				std::cout << "▉";

			Eigen::Vector3d vectorContinue(sampling.sample_points[ori]);  //与祖先Node的ori相同时容易出现问题，暂时先避免用相同ori
			int temp_now_last_node = now_last_node;
			bool jud_continue = false;

			for (int j = 0; j < save_ori[temp_now_last_node].size(); j++)
				if (vectorContinue == save_ori[temp_now_last_node][j]) {
					jud_continue = true;
					break;
				}

			if (jud_continue == true && ori != 0) {
				continue;
			}

			//rotating the blocks and then slicing//
			std::vector<Eigen::MatrixXd> temp_V;	//记录当前块的顶点信息，用于旋转
			temp_V.resize(Katana::Instance().temp_vertices.size());

			V_2.resize(V.rows());	//V_2记录输入模型网格V的顶点信息，用于旋转
			for (int i = 0; i < Katana::Instance().temp_vertices.size(); i++) {
				temp_V[i].resize(3, 1);
				temp_V[i](0, 0) = Katana::Instance().temp_vertices[i].x;
				temp_V[i](1, 0) = Katana::Instance().temp_vertices[i].y;
				temp_V[i](2, 0) = Katana::Instance().temp_vertices[i].z;
			}
			for (int i = 0; i < V.rows(); i++) {
				V_2[i].resize(3, 1);
				V_2[i](0, 0) = V.row(i).x();
				V_2[i](1, 0) = V.row(i).y();
				V_2[i](2, 0) = V.row(i).z();
			}
			Eigen::Vector3d vectorBefore(0, 0, 1);
			Eigen::Vector3d vectorAfter(sampling.sample_points[ori]);
			rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
			for (int i = 0; i < Katana::Instance().vertices.size(); i++)
				temp_V[i] = rotMatrix.inverse() * temp_V[i];	//将temp_V里的顶点坐标旋转，使得增材打印方向与z轴平行
			for (int i = 0; i < V.rows(); i++)
				V_2[i] = rotMatrix.inverse() * V_2[i];	//将V_2里的顶点坐标旋转，使得增材打印方向与z轴平行
			for (int i = 0; i < Katana::Instance().vertices.size(); i++) {
				Katana::Instance().vertices[i].x = temp_V[i](0, 0);
				Katana::Instance().vertices[i].y = temp_V[i](1, 0);
				Katana::Instance().vertices[i].z = temp_V[i](2, 0);	//将temp_V里的顶点坐标更新回Katana的vertices
			}

			for (int i = 0; i < Katana::Instance().triangles.size(); i++) {	//将每个三角形的顶点按z值从小到大排序，便于后续分层处理
				if (Katana::Instance().triangles[i].vertices[0]->z > Katana::Instance().triangles[i].vertices[1]->z) std::swap(Katana::Instance().triangles[i].vertices[0], Katana::Instance().triangles[i].vertices[1]);
				if (Katana::Instance().triangles[i].vertices[0]->z > Katana::Instance().triangles[i].vertices[2]->z) std::swap(Katana::Instance().triangles[i].vertices[0], Katana::Instance().triangles[i].vertices[2]);
				if (Katana::Instance().triangles[i].vertices[1]->z > Katana::Instance().triangles[i].vertices[2]->z) std::swap(Katana::Instance().triangles[i].vertices[1], Katana::Instance().triangles[i].vertices[2]);
			}
			start_time_6 = clock();
			vector<vector<vector<Vertex>>> all_slice_points;
			vector<vector<vector<Vertex>>> all_slice_points_contain;
			Katana::Instance().slicer.buildLayers();
			Katana::Instance().slicer.buildSegments();
			Katana::Instance().gcode.write(all_slice_points, all_slice_points_contain);	//katana分层并填充边界轮廓
			end_time_6 = clock();
			sum_time_4 += double(end_time_6 - start_time_6) / CLOCKS_PER_SEC;
			/////////////////////////////////////
			//std::cout << "aa" << Katana::Instance().vertices.size()<<" " << all_slice_points.size() << endl;
			//generate additive dependency graph//
			std::vector<Data> data;
			data.resize(1);
			data[0].ReadData(all_slice_points, all_slice_points_contain);
			Layer_Graph layer_graph(data[0]);
			//start_time_4 = clock();
			start_time_4 = clock();
			layer_graph.GetTrianglesForLayers(all_slice_points, Katana::Instance().map_segment_triangles, Katana::Instance().vertices, vectorAfter, height_of_beam_search, Tree_nodes_continue_id[now_last_node]);	//将切片轮廓映射到三角形集，建立每层三角形集合等中间信息？
			layer_graph.GenerateDependencyEdges();	//生成增材分层依赖图的边
			layer_graph.CollisionDetectionForAdditiveManufacturing(the_nozzle);	//增材的碰撞检测，标记增材的不可达点等信息
			end_time_4 = clock();
			sum_time_2 += double(end_time_4 - start_time_4) / CLOCKS_PER_SEC;

			//////////////////////////////////////
			vector<Eigen::MatrixXd> fragile_V;
			for (int i = 0; i < is_fragile_V_2[ori].size(); i++)
				if (is_fragile_V_2[ori][i] == true)
					fragile_V.push_back(V_2[i]);

			all_solutions_of_selected_layers.clear();
			all_solutions_of_selected_layers_contain.clear();
			all_cut_layers.clear(); all_cut_layers_dependency_layer.clear();
			pathes_include_S.clear(); pathes_include_sample_points.clear(); paths_include_covering_points.clear();
			//std::cout << "cc" << endl;

			bool flag_continue = false;	//是否继续从当前节点向下搜索
			bool previous_is_continue = false;	//记录上一个节点是否为continue节点
			/*if (Tree_nodes_judge_continue[now_last_node] == true)
				previous_is_continue = true;*/
			bool jud_admit = true;	//记录当前方向是否可行

			start_time_3 = clock();

			if (open_change_orientation == true)
				W1 = 1;
			DFS_search(layer_graph, flag_continue, previous_is_continue, judge_S_be_searched, judge_covering_points_be_searched, jud_admit);	//对当前方向进行增材分层依赖图的深度优先搜索，生成这个分层方向的解

			/*if (height_of_beam_search==1&& Tree_nodes_continue_id[now_last_node] >= 2)
				flag_continue = false;*/

			end_time_3 = clock();
			if (jud_admit == false)
				continue;
			sum_time += double(end_time_3 - start_time_3) / CLOCKS_PER_SEC;

			//cout << "()()()(" << double(end_time_3 - start_time_3) / CLOCKS_PER_SEC << endl;

			final_pathes_include_S.insert(final_pathes_include_S.end(), pathes_include_S.begin(), pathes_include_S.end());
			final_pathes_include_sample_points.insert(final_pathes_include_sample_points.end(), pathes_include_sample_points.begin(), pathes_include_sample_points.end());
			final_pathes_include_covering_points.insert(final_pathes_include_covering_points.end(), paths_include_covering_points.begin(), paths_include_covering_points.end());

			for (int i = 0; i < all_solutions_of_selected_layers.size(); i++) {	//将当前方向生成的所有解加入树结构中，作为下一层搜索的节点
				index_node++;
				Tree_nodes.push_back(all_solutions_of_selected_layers[i]);
				Tree_nodes_contain.push_back(all_solutions_of_selected_layers_contain[i]);
				Tree_nodes_cut_layers.push_back(all_cut_layers[i]);
				Tree_nodes_num_of_cut_layers_dependency_layer.push_back(all_cut_layers_dependency_layer[i]);
				Tree_nodes_fragile_V.push_back(fragile_V);
				//Tree_nodes_ori_all_the_area_S.push_back(ori_all_the_area_S);
				//Tree_nodes_judge_S_be_searched.push_back(judge_S_be_searched);
				pre_tree_nodes[index_node] = now_last_node;
				pre_index_of_nodes[index_node] = cont_number_of_queue;
				candidate_nodes.push_back(index_node);
				vector<Eigen::Vector3d> temp_vec;
				temp_vec.push_back(vectorAfter);
				save_ori.push_back(temp_vec);
				Tree_nodes_judge_continue.push_back(flag_continue);
				if (Tree_nodes_judge_continue[now_last_node] == true)
					Tree_nodes_continue_id.push_back(Tree_nodes_continue_id[now_last_node] + 1);
				else
					Tree_nodes_continue_id.push_back(0);
				Tree_nodes_error.push_back(false);
			}
		}
		end_time_2 = clock();


		if (Tree_nodes_judge_continue[now_last_node] == false)
			last_step_nodes.pop();
		cont_number_of_queue++;
		if (last_step_nodes.size() == 0 || Tree_nodes_judge_continue[now_last_node] == true) {
			cont_number_of_queue = 0;
			height_of_beam_search++;
			int cont_w = 0;  //应为0开始
			cout << endl << "Decomposing and sorting......" << endl;
			/////////sort_candidate_nodes//////////
			start_time = clock();
			vector<all_value> all_calculated_value;
			vector<all_value> pure_value;
			all_calculated_value.resize(candidate_nodes.size());
			pure_value.resize(candidate_nodes.size());
			for (int i = 0; i < candidate_nodes.size(); i++) {
				printf("\r[%d%%]>", i * 100 / (candidate_nodes.size() - 1));
				for (int j = 1; j <= i * 20 / candidate_nodes.size(); j++)
					cout << "▉";
				int index_of_pre_node = pre_index_of_nodes[candidate_nodes[i]];
				vector<vector<cv::Point3d>> all_cut_layers;
				vector<int> all_cut_layers_dependency_layer;
				all_cut_layers.clear(); all_cut_layers_dependency_layer.clear();
				for (int j = 0; j < Tree_nodes_cut_layers[candidate_nodes[i]].size(); j++) {
					int index_of_layers = Tree_nodes_cut_layers[candidate_nodes[i]][j];
					all_cut_layers.push_back(Tree_nodes[candidate_nodes[i]][index_of_layers]);
					all_cut_layers_dependency_layer.push_back(Tree_nodes_num_of_cut_layers_dependency_layer[candidate_nodes[i]][j]);
				}
				//clock_t start_time_8, end_time_8;

				Slicer_2 slicer_G;
				if (Tree_nodes_judge_continue[now_last_node] == false)
					slicer_G.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + ".obj");
				else
					slicer_G.load(file_name + "-" + to_string(height_of_beam_search - 1) + "_" + to_string(index_of_pre_node) + "_" + to_string(Tree_nodes_continue_id[candidate_nodes[i]] - 1) + "_subblock.obj");
				//start_time_8 = clock();
				all_calculated_value[i] = GainMesh(slicer_G, all_cut_layers, save_ori[candidate_nodes[i]][0], height_of_beam_search, cont_number_of_queue, index_of_pre_node, all_cut_layers_dependency_layer, Tree_nodes_judge_continue[now_last_node], Tree_nodes_continue_id[candidate_nodes[i]]);
				//end_time_8 = clock();
				//cout << "()()()(" << double(end_time_8 - start_time_8) / CLOCKS_PER_SEC << endl;
				//cout << endl <<candidate_nodes.size();

				calculate_fragile_value(all_calculated_value[i], all_cut_layers, Tree_nodes_fragile_V[candidate_nodes[i]]);

				Tree_nodes_larger_base.push_back(all_calculated_value[i].large_base);
				if (all_calculated_value[i].value_of_self_support == 0) {
					//cout << "is not self-support!" << endl;
					candidate_nodes.erase(candidate_nodes.begin() + i);
					all_calculated_value.erase(all_calculated_value.begin() + i);
					i--;
					continue;
				}

				//if(height_of_beam_search <=6)
					//cout << i << endl;

				detect_collision_with_printing_platform(i, candidate_nodes, all_calculated_value, all_cut_layers, save_ori[candidate_nodes[i]][0], the_nozzle);
			}

			sort_candidate_nodes(candidate_nodes, Tree_nodes, final_pathes_include_S, all_calculated_value, Tree_nodes_cut_layers, pre_tree_nodes, Tree_nodes_larger_base, final_pathes_include_covering_points, height_of_beam_search, save_ori, pure_value, Tree_nodes_continue_id[candidate_nodes[0]]);
			end_time = clock();

			ofile_cont << final_pathes_include_S[candidate_nodes[0]].size() << endl;
			double temp_record_size = 0;
			for (int j = 0; j < Tree_nodes[candidate_nodes[0]].size(); j++) {
				for (int k = 0; k < Tree_nodes[candidate_nodes[0]][j].size() - 1; k++) {
					temp_record_size += distance3d(Tree_nodes[candidate_nodes[0]][j][k], Tree_nodes[candidate_nodes[0]][j][k + 1]);
				}
			}
			ofile_cont_size << temp_record_size << endl;
			//////////////////////////////////////
			if (Tree_nodes_judge_continue[now_last_node] == true) {
				W1 = 1;
			}

			//cout << "bbb";

			bool jud_continue_last_node = false;
			//	for (int i = 0; i < candidate_nodes.size(); i++)
			//		cout << save_ori[candidate_nodes[i]][0].x() << " " << save_ori[candidate_nodes[i]][0].y() << " " << save_ori[candidate_nodes[i]][0].z() << " " << endl;
			start_time_7 = clock();
			while (candidate_nodes.size() != 0 && cont_w < W1 && cont_w < candidate_nodes.size()) {
				cout << endl << "selected orientation: " << save_ori[candidate_nodes[cont_w]][0].x() << " " << save_ori[candidate_nodes[cont_w]][0].y() << " " << save_ori[candidate_nodes[cont_w]][0].z();
				cout << endl << "number of candidate_nodes: " << candidate_nodes.size() << endl;
				Sum_candidate_blocks += candidate_nodes.size();
				Sum_connected_components += Tree_nodes_cut_layers[candidate_nodes[cont_w]].size();
				if (Tree_nodes_judge_continue[now_last_node] == false)
					last_step_nodes.push(candidate_nodes[cont_w]);
				cout << "self support value of selected node: " << all_calculated_value[cont_w].value_of_self_support << endl << "■■■■■■■■■■■■■■■■■■■" << endl;
				cout << "value_of_more_slice_layers: " << all_calculated_value[cont_w].value_of_more_slice_layers << endl;
				cout << "value_of_area_S: " << all_calculated_value[cont_w].value_of_area_S << endl;
				cout << "value_of_covering_points: " << all_calculated_value[cont_w].value_of_covering_points << endl;
				cout << "value_of_less_clipping_plane: " << all_calculated_value[cont_w].value_of_less_clipping_plane << endl;
				cout << "value_of_fragile: " << all_calculated_value[cont_w].value_of_fragile << endl;
				cout << "value_of_orientation: " << all_calculated_value[cont_w].value_of_orientation << endl;
				cout << "value_of_projected: " << all_calculated_value[cont_w].value_of_projected << endl;
				cout << "■■■■■■■■■■■■■■■■■■■" << endl;

				cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
				cout << "Pure value orien:" << pure_value[cont_w].value_of_orientation << endl;
				cout << "Pure value fragile:" << pure_value[cont_w].value_of_fragile << endl;
				cout << "Pure value projection:" << pure_value[cont_w].value_of_projected << endl;
				cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
				if (Tree_nodes_judge_continue[candidate_nodes[cont_w]] == true) {
					cout << "*****NEED CHANGE ORIENTATION*****" << endl;
					cont_extra_additive_orientation++;
				}

				evaluation_value[0] += all_calculated_value[cont_w].value_of_more_slice_layers;
				evaluation_value[1] += all_calculated_value[cont_w].value_of_covering_points;
				evaluation_value[2] += all_calculated_value[cont_w].value_of_less_clipping_plane;
				evaluation_value[3] += all_calculated_value[cont_w].value_of_orientation;
				evaluation_value[4] += all_calculated_value[cont_w].value_of_fragile;
				evaluation_value[5] += all_calculated_value[cont_w].value_of_projected;
				//decompose the model for every node//
				int index_of_pre_node = pre_index_of_nodes[candidate_nodes[cont_w]];
				vector<vector<cv::Point3d>> all_cut_layers;
				vector<int> flag_cut_layers_is_hole;
				vector<int> all_cut_layers_dependency_layer;
				all_cut_layers.clear(); all_cut_layers_dependency_layer.clear();
				flag_cut_layers_is_hole.clear();
				for (int j = 0; j < Tree_nodes_cut_layers[candidate_nodes[cont_w]].size(); j++) {
					int index_of_layers = Tree_nodes_cut_layers[candidate_nodes[cont_w]][j];
					all_cut_layers.push_back(Tree_nodes[candidate_nodes[cont_w]][index_of_layers]);
					flag_cut_layers_is_hole.push_back(-1);
					if (Tree_nodes_contain[candidate_nodes[cont_w]][index_of_layers].size() != 0) {
						all_cut_layers.push_back(Tree_nodes_contain[candidate_nodes[cont_w]][index_of_layers]);
						all_cut_layers_dependency_layer.push_back(Tree_nodes_num_of_cut_layers_dependency_layer[candidate_nodes[cont_w]][j]);
						flag_cut_layers_is_hole.push_back(all_cut_layers.size() - 2);
					}
					all_cut_layers_dependency_layer.push_back(Tree_nodes_num_of_cut_layers_dependency_layer[candidate_nodes[cont_w]][j]);
				}

				cout << "***** number_of_S: " << final_pathes_include_S[candidate_nodes[cont_w]].size() << endl;
				cout << "***** all_cut_layers.size(): " << all_cut_layers.size() << endl;

				std::vector<TRiangle> current_remove_triangles;
				Slicer_2 current_slicer;
				//cout << "aaa";
				/*for (int i = 0; i < Tree_nodes_continue_id.size(); i++)
					cout << "*"<<Tree_nodes_continue_id[i] << endl;*/

				CutMesh(Tree_nodes[candidate_nodes[cont_w]], Tree_nodes_contain[candidate_nodes[cont_w]],
					all_cut_layers, save_ori[candidate_nodes[cont_w]][0],
					height_of_beam_search,
					cont_number_of_queue,
					index_of_pre_node,
					all_cut_layers_dependency_layer,
					jud_outer_beam_search_terminate,
					current_remove_triangles,
					current_slicer,
					Tree_nodes_judge_continue[candidate_nodes[cont_w]], Tree_nodes_judge_continue[now_last_node],
					pre_index_of_nodes[now_last_node],
					Tree_nodes_error,
					candidate_nodes[cont_w],
					Tree_nodes_continue_id[candidate_nodes[cont_w]],
					flag_cut_layers_is_hole);


				//if (Tree_nodes_judge_continue[now_last_node] == true) {
				//	cout << "aaa";
				//	Tree_nodes.pop_back();
				//	Tree_nodes_cut_layers.pop_back();
				//	Tree_nodes_num_of_cut_layers_dependency_layer.pop_back();
				//	Tree_nodes_fragile_V.pop_back();
				//	index_node--;
				//	candidate_nodes.pop_back();
				//	save_ori.pop_back();
				//	Tree_nodes_judge_continue.pop_back();
				//	Tree_nodes_error.pop_back();
				//	height_of_beam_search--;
				//	Tree_nodes_judge_continue[now_last_node] = false;
				//	//Tree_nodes_continue_id[now_last_node]++;
				//}
				cout << "&&&&" << height_of_beam_search << endl;
				if (Tree_nodes_judge_continue[candidate_nodes[0]] == true) {
					cout << "AAAA!" << endl;
					jud_continue_last_node = true;

				}
				if (Tree_nodes_judge_continue[now_last_node] == true) { //Tree_nodes_judge_continue[now_last_node] == true  //应该改为Tree_nodes_judge_continue[candidate_nodes[0]] == true?
					cout << "BBBB!" << candidate_nodes[0] << endl;
					if (Tree_nodes_judge_continue[candidate_nodes[0]] == false)
						Tree_nodes_judge_continue[now_last_node] = false;
					height_of_beam_search--;
					Tree_nodes.pop_back();
					Tree_nodes_contain.pop_back();
					Tree_nodes_cut_layers.pop_back();
					Tree_nodes_num_of_cut_layers_dependency_layer.pop_back();
					Tree_nodes_fragile_V.pop_back();
					index_node--;
					candidate_nodes.pop_back();
					Tree_nodes_judge_continue.pop_back();
					Tree_nodes_error.pop_back();
					//Tree_nodes_judge_continue[now_last_node] = false;
					Tree_nodes_continue_id[now_last_node] = Tree_nodes_continue_id[candidate_nodes[0]];
					save_ori[now_last_node].push_back(save_ori[candidate_nodes[0]][0]);
					Tree_nodes_continue_id.pop_back();
					save_ori.pop_back();

				}


				if (height_of_beam_search <= 15) {
					//cout << save_ori[candidate_nodes[cont_w]].x() <<" "<< save_ori[candidate_nodes[cont_w]].y() << " " << save_ori[candidate_nodes[cont_w]].z() << endl;
					subtractive_accessibility_decomposition(current_remove_triangles, height_of_beam_search, index_of_pre_node, cutting_tool, current_slicer);
				}

				subtractive_remove_output(current_remove_triangles, current_slicer, height_of_beam_search);

				cont_number_of_queue++;
				//////////////////////////////////////

				cont_w++;
			}
			end_time_7 = clock();
			std::cout << "&&&time&&& Evaluation " << double(end_time - start_time + end_time_7 - start_time_7) / CLOCKS_PER_SEC << std::endl;
			std::cout << "&&&time&&& Update subtractive dependency graph: " << sum_time_3 << std::endl;
			std::cout << "&&&time&&& Slicing: " << sum_time_4 << std::endl;
			std::cout << "&&&time&&& Build addictive dependency graph: " << sum_time_2 << std::endl;
			std::cout << "&&&time&&& Co graph merging: " << sum_time << std::endl;
			total_time_1 += sum_time_3;
			total_time_2 += sum_time_2 + sum_time_4 + sum_time_5;
			total_time_3 += sum_time;
			total_time_4 += double(end_time - start_time + end_time_7 - start_time_7) / CLOCKS_PER_SEC;
			sum_time = 0;
			sum_time_2 = 0;
			sum_time_3 = 0;
			sum_time_4 = 0;
			sum_time_5 = 0;
			//std::cout << "&&&time&&& a step of beam search: " << double(end_time_2 - start_time_2) / CLOCKS_PER_SEC << std::endl;
			/*if (height_of_beam_search == 5)
				W1 = 20;*/
				/*if (height_of_beam_search >= 2 && Tree_nodes_continue_id[candidate_nodes[0]] >= 0) {
					end_time_total = clock();
					std::cout << "&&&time&&& total of outer-beam search: " << double(end_time_total - start_time_total) / CLOCKS_PER_SEC << "s &&&&&&&" << std::endl << endl;
					system("pause");
				} */

				//system("pause");
				//save_ori.clear();
			cont_number_of_queue = 0;
			if (jud_outer_beam_search_terminate == true) {
				/*cout << "$$$$$$$$$$$$$$$$$$$$$$$$$ Data $$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
				cout << "#S: " << all_voronoi_cells.size() << endl;
				cout << "#P_i: " << num_inaccessible_points << endl;
				cout << "#R_p: " << float(num_inaccessible_points) / float(all_voronoi_cells.size()) << endl;
				cout << "#B_C !!: " << Sum_candidate_blocks << endl;
				cout << "#I: " << height_of_beam_search << endl;
				cout << "#B: " << height_of_beam_search << endl;
				cout << "#O: " << cont_extra_additive_orientation + height_of_beam_search << endl;
				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$ Total time $$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
				cout << "Build subtractive graph: " << time_build_subtractive_graph << endl;
				cout << "Update subtractive dependency graph: " << total_time_1 << endl;
				cout << "Build addictive dependency graph: " << total_time_2 << endl;
				cout << "Co graph merging: " << total_time_3 << endl;
				cout << "Evaluation: " << total_time_4 << endl;
				cout << "Total: " << time_build_subtractive_graph + total_time_1 + total_time_2 + total_time_3 + total_time_4 << endl;
				cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;*/

				double global_score = 1 / (height_of_beam_search);
				fstream f;
				//追加写入,在原来基础上加了ios::app  faith
				std::string txt_name = file_name + ".txt";
				f.open(txt_name.c_str(), ios::out | ios::app);
				//输入你想写入的内容 
				f << global_score << endl;
				f.close();

				break;
			}
			if (height_of_beam_search >= 10) {
				fstream f;
				//追加写入,在原来基础上加了ios::app  faith
				std::string txt_name = file_name + ".txt";
				f.open(txt_name.c_str(), ios::out | ios::app);
				//输入你想写入的内容 
				f << 0 << endl;
				f.close();
				break;
			}
			if (jud_continue_last_node == true)
				cout << endl << "--------------------- height of outer-beam search: " << height_of_beam_search << " ---------------------";
			else
				cout << endl << "--------------------- height of outer-beam search: " << height_of_beam_search + 1 << " ---------------------";
			candidate_nodes.clear();
		}
	}

	vector<int> last_available_block;
	vector<int> exist_points;
	vector<vector<int>> is_point_exist_in_block;
	vector<vector<int>> is_point_available_in_block;
	vector<int> is_available(height_of_beam_search + 1, 1);


	last_available_block.resize(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		last_available_block[i] = height_of_beam_search;
	}

	for (int i = 0; i < V.rows(); i++) {
		is_point_exist_in_block.push_back(is_available);
	}

	for (int i = 0; i < V.rows(); i++)
		is_point_available_in_block.push_back(is_available);

	cout << "TTYTYTYTYTY";

	//subtractive_accessibility_decomposition_within_2_blocks(height_of_beam_search, cutting_tool);
}

void HybridManufacturing::DFS_search(Layer_Graph layer_graph, bool& flag_continue, bool previous_is_continue, vector<bool> judge_S_be_searched, vector<bool> judge_covering_points_be_searched, bool& jud_admit)
{
	int W2 = 1;  //2
	std::vector<int> Tree_nodes;
	vector<vector<int>> Tree_nodes_for_S;
	vector<vector<int>> Tree_nodes_for_sample_points;
	vector<vector<int>> Tree_nodes_for_covering_points;
	int* pre_tree_nodes;
	pre_tree_nodes = new int[100000];
	vector<int> candidate_nodes;
	queue<int> last_step_nodes;
	vector<int> S_in_block;
	vector<int> terminate_nodes;
	memset(pre_tree_nodes, -1, sizeof(pre_tree_nodes));
	int root_node = -1;
	Tree_nodes.push_back(root_node);
	vector<int> root_node_for_S;
	Tree_nodes_for_S.push_back(root_node_for_S);
	Tree_nodes_for_sample_points.push_back(root_node_for_S);
	Tree_nodes_for_covering_points.push_back(root_node_for_S);
	int index_node = 0;
	last_step_nodes.push(index_node);
	//vector<vector<area_S>> ori_all_the_area_S = all_the_area_S;

	//std::cout << "AA" << endl;
	//find boundary point of every layer//
	vector<double> left_point_of_layer, right_point_of_layer, top_point_of_layer, bottom_point_of_layer;
	left_point_of_layer.resize(layer_graph.total_node_num); right_point_of_layer.resize(layer_graph.total_node_num);
	top_point_of_layer.resize(layer_graph.total_node_num); bottom_point_of_layer.resize(layer_graph.total_node_num);

	vector<vector<Point_2>> vec_points_2d(layer_graph.total_node_num);
	vector<Polygon_2> vec_polygon(layer_graph.total_node_num);
	/*Point_2 pp(V(i, 0), V(i, 1));
	if (polygon.bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE)
		jud_still_exist = true;*/

	for (int i = 0; i < layer_graph.total_node_num; i++) {
		left_point_of_layer[i] = MAX_I;
		right_point_of_layer[i] = -MAX_I;
		top_point_of_layer[i] = -MAX_I;
		bottom_point_of_layer[i] = MAX_I;
	}
	int temp_num = 0;
	for (int i = 0; i < layer_graph.data.slice_points.size(); i++) {
		for (int j = 0; j < layer_graph.data.slice_points[i].size(); j++) {
			for (int k = 0; k < layer_graph.data.slice_points[i][j].size(); k++) {
				/*if (layer_graph.data.slice_points[i][j][k].x < left_point_of_layer[temp_num])
					left_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].x;
				if (layer_graph.data.slice_points[i][j][k].x > right_point_of_layer[temp_num])
					right_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].x;*/
					/*if (layer_graph.data.slice_points[i][j][k].y < bottom_point_of_layer[temp_num])
						bottom_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].y;
					if (layer_graph.data.slice_points[i][j][k].y > top_point_of_layer[temp_num])
						top_point_of_layer[temp_num] = layer_graph.data.slice_points[i][j][k].y;*/
             vec_points_2d[temp_num].push_back(Point_2(layer_graph.data.slice_points[i][j][k].x(), layer_graph.data.slice_points[i][j][k].y()));
			}
			vec_polygon[temp_num] = Polygon_2(vec_points_2d[temp_num].begin(), vec_points_2d[temp_num].end());
			//std::cout << "!!@@##" << std::endl;
			if (!vec_polygon[temp_num].is_simple()) {

				vec_polygon[temp_num] = Polygon_2(vec_points_2d[temp_num].begin(), vec_points_2d[temp_num].end() - 1);

				if (!vec_polygon[temp_num].is_simple()) {
					std::cout << "polygon is not simple HybridManufacturing::DFS_search cutcut!" << std::endl;
				}
			}
			temp_num++;
		}
	}


	//////////////////////////////////////
	//std::cout << "BB" << endl;
	//find all the sample points in every layer//
	vector<vector<int>> sample_point_in_layer;	//存储每个layer中包含的sample point的索引
	sample_point_in_layer.resize(layer_graph.data.total_node_num);
	vector<bool> judge_sample_point_be_searched;
	judge_sample_point_be_searched.resize(V_2.size());
	for (int i = 0; i < V_2.size(); i++)
		judge_sample_point_be_searched[i] = false;
	for (int i = 0; i < layer_graph.data.total_node_num; i++) {
		for (int j = 0; j < V_2.size(); j++) {
			pair<int, int> index_slice_layer = layer_graph.data.index[i];
			Point_2 pp(V_2[j](0, 0), V_2[j](1, 0));
			if (judge_sample_point_be_searched[j] == false
				&& vec_polygon[i].bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE // bounded_sideh很耗时?
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[j](2, 0) > 0
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[j](2, 0) <= dh) {
				sample_point_in_layer[i].push_back(j);	//V_2的第j个点在第i个layer中
				judge_sample_point_be_searched[j] = true;	//标记该点需要被搜索
			}
		}
	}
	////////////////////////////////////////

	//std::cout << "CC" << endl;
	//find area S and covering points in the node//
	vector<vector<int>> temp_all_S_in_the_block, temp_all_covering_points_in_the_block;
	temp_all_S_in_the_block.clear(); temp_all_covering_points_in_the_block.clear();
	temp_all_S_in_the_block.resize(layer_graph.data.total_node_num);
	temp_all_covering_points_in_the_block.resize(layer_graph.data.total_node_num);
	for (int i = 0; i < layer_graph.data.total_node_num; i++) {
		pair<int, int> index_slice_layer = layer_graph.data.index[i];
		for (int j = 0; j < ori_num_points_of_ori_in_all_the_area_S.size(); j++) {
			Point_2 pp(V_2[map_S_and_vertex[j]](0, 0), V_2[map_S_and_vertex[j]](1, 0));
			if (judge_S_be_searched[j] == false
				&& vec_polygon[i].bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[map_S_and_vertex[j]](2, 0) > 0
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[map_S_and_vertex[j]](2, 0) <= dh) {
				temp_all_S_in_the_block[i].push_back(j);
				judge_S_be_searched[j] = true;	//标记该area S需要被搜索
			}
		}
		for (int j = 0; j < ori_all_the_covering_points.size(); j++) {
			Point_2 pp(V_2[map_covering_points_and_vertex[j]](0, 0), V_2[map_covering_points_and_vertex[j]](1, 0));
			if (judge_covering_points_be_searched[j] == false
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[map_covering_points_and_vertex[j]](2, 0) > 0
				&& layer_graph.data.z_value[index_slice_layer.first][index_slice_layer.second][0] - V_2[map_covering_points_and_vertex[j]](2, 0) <= dh
				&& vec_polygon[i].bounded_side(pp) != CGAL::ON_UNBOUNDED_SIDE) {
				temp_all_covering_points_in_the_block[i].push_back(j);
				judge_covering_points_be_searched[j] = true;	//标记该covering point需要被搜索
			}
		}
	}
	////////////////////////////////////////
	//std::cout << "DD" << endl;
	while (last_step_nodes.size() != 0 || candidate_nodes.size() != 0) {
		flag_continue = false;
		int now_last_node = last_step_nodes.front();
		bool jud_terminate = true;

		//update//
		while (pre_tree_nodes[now_last_node] != -1) {
			layer_graph.UpdateDegree_2(Tree_nodes[now_last_node], -1);
			layer_graph.node_visited[Tree_nodes[now_last_node]] = true;
			now_last_node = pre_tree_nodes[now_last_node];
		}
		//////////
		now_last_node = last_step_nodes.front();

		for (int i = 0; i < layer_graph.total_node_num; i++) {
			int v = i;
			bool flag_self_support = true;
			if (Tree_nodes.size() != 0 && Tree_nodes[now_last_node] == v) continue;
			bool jud_continue_2 = false;
			/*for (int j = 0; j < candidate_nodes.size(); j++)
				if (Tree_nodes[candidate_nodes[j]] == v) {
					jud_continue_2 = true;
					break;
				}*/
			if (jud_continue_2 == true) continue;
			if (layer_graph.node_visited[v]) continue;
			if (layer_graph.out_degree[v] != 0) continue;   //out-degree == 0
			//////////collision dependency edges///////////////
			if (Tree_nodes.size() != 0 && layer_graph.IsDepend_collision(Tree_nodes[now_last_node], v))
				continue;
			///////////////////////////////////////////////////
			if (layer_graph.is_the_layer_self_suppot[v] == false && now_last_node != 0) {  //self-support constraint
				flag_self_support = false;
			}

			vector<int> all_S_in_the_block, all_covering_points_in_the_block;
			all_S_in_the_block = temp_all_S_in_the_block[i];
			all_covering_points_in_the_block = temp_all_covering_points_in_the_block[i];

			//try merge//
			bool jud_merge = true;
			for (int j = 0; j < all_S_in_the_block.size(); j++) {
				//check whether all S in the block is accessible//
				bool jud_merge_2 = false;
				for (int k = 0; k < ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]].size(); k++) {
					if (ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]][k] < 0)
						cout << "error!!";
					if (ori_num_points_of_ori_in_all_the_area_S[all_S_in_the_block[j]][k] == 0) {
						jud_merge_2 = true;
						break;
					}
				}
				if (jud_merge_2 == false)
					jud_merge = false;
				//////////////////////////////////////////////////
			}
			if (jud_merge == true && flag_self_support == false && previous_is_continue == false) {
				if (open_change_orientation == true)
					flag_continue = true;
			}
			if (jud_merge == true && flag_self_support == true) {
				jud_terminate = false;
				index_node++;
				Tree_nodes.push_back(v);

				Tree_nodes_for_S.push_back(all_S_in_the_block);
				Tree_nodes_for_sample_points.push_back(sample_point_in_layer[v]);
				Tree_nodes_for_covering_points.push_back(all_covering_points_in_the_block);
				pre_tree_nodes[index_node] = now_last_node;
				candidate_nodes.push_back(index_node);

				break; //DFS
			}
			////////////
		}
		if (jud_terminate == true) {
			terminate_nodes.push_back(now_last_node);
			break;
		}

		//update//
		while (pre_tree_nodes[now_last_node] != -1) {
			layer_graph.UpdateDegree_2(Tree_nodes[now_last_node], 1);
			layer_graph.node_visited[Tree_nodes[now_last_node]] = false;
			now_last_node = pre_tree_nodes[now_last_node];
		}
		//////////

		//sort_candidate_nodes(candidate_nodes, Tree_nodes_for_S);
		last_step_nodes.pop();
		//cout << candidate_nodes.size() << endl;
		if (last_step_nodes.size() == 0) {
			int cont_w = 0;
			while (candidate_nodes.size() != 0 && cont_w < W2 && cont_w < candidate_nodes.size()) {
				last_step_nodes.push(candidate_nodes[cont_w]);
				cont_w++;
			}
			candidate_nodes.clear();
		}
	}

	//every path (blcok) as a candidate node of outer-beam search (save slice layers)//
	//final_pathes save the index of slice layer (layer_graph.total_node_num)
	pathes_include_S.clear();
	pathes_include_sample_points.clear();
	paths_include_covering_points.clear();
	vector<vector<int>> final_pathes;
	vector<int> final_pathes_height;
	for (int i = 0; i < terminate_nodes.size(); i++) {
		int height = 0;
		vector<int> current_path, current_path_include_S, current_path_include_sample_points, current_path_include_covering_points;
		current_path.clear();
		current_path_include_S.clear();
		current_path_include_sample_points.clear();
		current_path_include_covering_points.clear();
		int current_node = terminate_nodes[i];
		current_path.push_back(Tree_nodes[terminate_nodes[i]]);
		for (int j = 0; j < Tree_nodes_for_S[terminate_nodes[i]].size(); j++)
			current_path_include_S.push_back(Tree_nodes_for_S[terminate_nodes[i]][j]);
		for (int j = 0; j < Tree_nodes_for_sample_points[terminate_nodes[i]].size(); j++)
			current_path_include_sample_points.push_back(Tree_nodes_for_sample_points[terminate_nodes[i]][j]);
		for (int j = 0; j < Tree_nodes_for_covering_points[terminate_nodes[i]].size(); j++)
			current_path_include_covering_points.push_back(Tree_nodes_for_covering_points[terminate_nodes[i]][j]);
		while (pre_tree_nodes[current_node] > 0) {
			current_node = pre_tree_nodes[current_node];
			current_path.push_back(Tree_nodes[current_node]);
			height++;
			for (int j = 0; j < Tree_nodes_for_S[current_node].size(); j++) {
				current_path_include_S.push_back(Tree_nodes_for_S[current_node][j]);
			}

			for (int j = 0; j < Tree_nodes_for_sample_points[current_node].size(); j++)
				current_path_include_sample_points.push_back(Tree_nodes_for_sample_points[current_node][j]);
			for (int j = 0; j < Tree_nodes_for_covering_points[current_node].size(); j++)
				current_path_include_covering_points.push_back(Tree_nodes_for_covering_points[current_node][j]);
		}
		std::reverse(current_path.begin(), current_path.end());
		std::reverse(current_path_include_S.begin(), current_path_include_S.end());
		std::reverse(current_path_include_sample_points.begin(), current_path_include_sample_points.end());
		std::reverse(current_path_include_covering_points.begin(), current_path_include_covering_points.end());
		final_pathes.push_back(current_path);
		final_pathes_height.push_back(height);
		pathes_include_S.push_back(current_path_include_S);
		pathes_include_sample_points.push_back(current_path_include_sample_points);
		paths_include_covering_points.push_back(current_path_include_covering_points);
	}

	//sort final_pathes and select W2 pathes//
	for (int i = 0; i < final_pathes.size(); i++)
		for (int j = i + 1; j < final_pathes.size(); j++) {
			if (final_pathes_height[i] < final_pathes_height[j]) {
				swap(final_pathes[i], final_pathes[j]);
				swap(final_pathes_height[i], final_pathes_height[j]);
			}
		}
	while (final_pathes.size() > W2) {
		auto itr = final_pathes.begin() + W2;
		final_pathes.erase(itr);
		auto itr2 = pathes_include_S.begin() + W2;
		pathes_include_S.erase(itr2);
		auto itr3 = pathes_include_sample_points.begin() + W2;
		pathes_include_sample_points.erase(itr3);
		auto itr4 = paths_include_covering_points.begin() + W2;
		paths_include_covering_points.erase(itr4);
	}
	////////////////////

	all_cut_layers = FindAllCutLayers(layer_graph, final_pathes, all_cut_layers_dependency_layer, jud_admit);

	vector<vector<vector<cv::Point3d>>> all_selected_points;	//存储所有被选择的路径上的切片轮廓的各个点的坐标
	vector<vector<vector<cv::Point3d>>> all_selected_points_contain;
	for (int i = 0; i < final_pathes.size(); i++) {
		vector<vector<cv::Point3d>> temp_vec_1;
		all_selected_points.push_back(temp_vec_1);
		all_selected_points_contain.push_back(temp_vec_1);
		for (int j = 0; j < final_pathes[i].size(); j++) {
     vector<cv::Point3d> temp_vec_2;
			all_selected_points[i].push_back(temp_vec_2);
			all_selected_points_contain[i].push_back(temp_vec_2);
			pair<int, int> index_slice_point = layer_graph.data.index[final_pathes[i][j]];
			for (int k = 0; k < layer_graph.data.slice_points[index_slice_point.first][index_slice_point.second].size(); k++) {
                cv::Point3d current_point(layer_graph.data.slice_points[index_slice_point.first][index_slice_point.second][k].x(),
					layer_graph.data.slice_points[index_slice_point.first][index_slice_point.second][k].y(),
					layer_graph.data.z_value[index_slice_point.first][index_slice_point.second][k]);
				all_selected_points[i][j].push_back(current_point);
			}
			for (int k = 0; k < layer_graph.data.slice_points_contain[index_slice_point.first][index_slice_point.second].size(); k++) {
                cv::Point3d current_point(layer_graph.data.slice_points_contain[index_slice_point.first][index_slice_point.second][k].x(),
					layer_graph.data.slice_points_contain[index_slice_point.first][index_slice_point.second][k].y(),
					layer_graph.data.z_value[index_slice_point.first][index_slice_point.second][0]);
				all_selected_points_contain[i][j].push_back(current_point);
			}
		}
	}


	//////////visualization all triangles cross by layers//////////
	//vector<vector<vector<Vertex*>>> show_triangles(final_pathes[0].size());
	//for (int i = 0; i < 1; i++) {
	//	for (int j = 0; j < final_pathes[0].size(); j++) {
	//		show_triangles[j].resize(layer_graph.all_triangles_of_layers[final_pathes[i][j]].size());
	//		for (int k = 0; k < layer_graph.all_triangles_of_layers[final_pathes[i][j]].size(); k++) {
	//			Vertex* v1 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[0];
	//			Vertex* v2 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[1];
	//			Vertex* v3 = layer_graph.all_triangles_of_layers[final_pathes[i][j]][k]->vertices_2[2];
	//			//double ans = (v2->x - v1->x) * (v2->y - v3->y) - (v2->y - v1->y) * (v2->x - v3->x);
	//			//if (ans > 0)	//is clockwise
	//			//	swap(v2, v3);
	//			show_triangles[j][k].push_back(v1);
	//			show_triangles[j][k].push_back(v2);
	//			show_triangles[j][k].push_back(v3);
	//		}
	//	}
	//}
	//std::ofstream dstream(".\\vis\\show_triangles.stl");
	//if (!dstream.is_open()) {
	//	std::cout << "can not open " << std::endl;
	//}
	//dstream << "solid STL generated by MeshLab" << std::endl;
	//for (int j = 0; j < final_pathes[0].size(); j++) {
	//	for (int k = 0; k < show_triangles[j].size(); k++) {
	//		dstream << "  facet normal " << "0 0 0" << std::endl;
	//		dstream << "    outer loop" << std::endl;
	//		for (int l = 0; l < 3; l++) {
	//			dstream << "      vertex  " << show_triangles[j][k][l]->x << " " << show_triangles[j][k][l]->y << " " << show_triangles[j][k][l]->z << std::endl;
	//		}
	//		dstream << "    endloop" << std::endl;
	//		dstream << "  endfacet" << std::endl;
	//	}
	//}
	//dstream << "endsolid vcg" << std::endl;
	//dstream.close();

	delete[] pre_tree_nodes;

	all_solutions_of_selected_layers = all_selected_points;
	all_solutions_of_selected_layers_contain = all_selected_points_contain;
	return;
	//////////////////////////////////////////////////////////////////////////////////
}

void HybridManufacturing::sort_candidate_nodes(vector<int>& candidate_nodes, vector<vector<vector<cv::Point3d>>> Tree_nodes, vector<vector<int>> final_pathes_include_S, vector<all_value>& all_calculated_value, vector<vector<int>> Tree_nodes_cut_layers, int pre_tree_nodes[], vector<double> Tree_nodes_larger_base, vector<vector<int>> final_pathes_include_covering_points, int height_of_beam_search, vector<vector<Eigen::Vector3d>> save_ori, vector<all_value>& pure_value, int id_continue)
{
	double W_less_area_S = 0, W_more_slice_layers = 0.6, W_covering_points = 0.4, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;

	//double W_less_area_S = 0, W_more_slice_layers = 0.1, W_covering_points = 0.1, W_less_clipping_plane = 0.4, W_fragile = 0, W_larger_base = 0, W_orientation = 0.4, W_projected = 0;
	//double W_less_area_S = 0, W_more_slice_layers = 0, W_covering_points = 0, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 1;

	//if(height_of_beam_search >= 12)  //test for figure of paper
	//	W_less_area_S = 0.1, W_more_slice_layers = 1, W_covering_points = 0.1, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;
	/*if(height_of_beam_search == 2)
		W_less_area_S = 0.9, W_more_slice_layers = 0.1, W_covering_points = 0, W_less_clipping_plane = 0, W_fragile = 0, W_larger_base = 0, W_orientation = 0, W_projected = 0;*/

		//random test
		//int rand_id = rand() % all_calculated_value.size();
		/*int rand_id = 5;
		swap(candidate_nodes[0], candidate_nodes[rand_id]);
		swap(all_calculated_value[0], all_calculated_value[rand_id]);
		return;*/

	for (int i = 0; i < all_calculated_value.size(); i++) {
		if (all_calculated_value[i].number_of_remaining_face < terminate_threshold_of_number_of_faces && save_ori[candidate_nodes[i]][0].x() == 0 && save_ori[candidate_nodes[i]][0].y() == 0) {
			swap(candidate_nodes[0], candidate_nodes[i]);
			swap(all_calculated_value[0], all_calculated_value[i]);
			cout << "Remaining face: " << all_calculated_value[i].number_of_remaining_face << endl;
			return;
		}
	}
	for (int i = 0; i < all_calculated_value.size(); i++) {
		if (height_of_beam_search == 4 && save_ori[candidate_nodes[i]][0].x() == 0 && save_ori[candidate_nodes[i]][0].y() == 0) {
			swap(candidate_nodes[0], candidate_nodes[i]);
			swap(all_calculated_value[0], all_calculated_value[i]);
			cout << "Remaining face: " << all_calculated_value[i].number_of_remaining_face << endl;
			return;
		}
	}

	/*if(height_of_beam_search == 2)
		W_less_area_S = 0, W_more_slice_layers = 0.1, W_covering_points = 0.1, W_less_clipping_plane = 0.8, W_fragile = 0, W_larger_base = 0;*/

		//normalization - local//
		/*int max_area_S = -MAX_I, min_area_S = MAX_I;
		int max_slice_layers = -MAX_I, min_slice_layers = MAX_I;
		int max_clipping_plane = -MAX_I, min_clipping_plane = MAX_I;
		double max_large_base = -MAX_I, min_large_base = MAX_I;
		for (int i = 0; i < all_calculated_value.size(); i++) {
			max_large_base = max(max_large_base, all_calculated_value[i].large_base);
			min_large_base = min(min_large_base, all_calculated_value[i].large_base);

			max_area_S = max(max_area_S, int(final_pathes_include_S[candidate_nodes[i]].size()));
			min_area_S = min(min_area_S, int(final_pathes_include_S[candidate_nodes[i]].size()));

			max_slice_layers = max(max_slice_layers, int(Tree_nodes[candidate_nodes[i]].size()));
			min_slice_layers = min(min_slice_layers, int(Tree_nodes[candidate_nodes[i]].size()));

			max_clipping_plane = max(max_clipping_plane, int(Tree_nodes_cut_layers[candidate_nodes[i]].si
			ze()));
		}


		for (int i = 0; i < all_calculated_value.size(); i++) {
			if (max_large_base - min_large_base != 0)
				all_calculated_value[i].large_base = (all_calculated_value[i].large_base - min_large_base) / (max_large_base - min_large_base);
			else
				all_calculated_value[i].large_base = 0;

			if (max_area_S - min_area_S != 0)
				all_calculated_value[i].value_of_area_S = double(final_pathes_include_S[candidate_nodes[i]].size() - min_area_S) / double(max_area_S - min_area_S);
			else
				all_calculated_value[i].value_of_area_S = 0;

			if (max_slice_layers - min_slice_layers != 0)
				all_calculated_value[i].value_of_more_slice_layers = double(Tree_nodes[candidate_nodes[i]].size() - min_slice_layers) / double(max_slice_layers - min_slice_layers);
			else
				all_calculated_value[i].value_of_more_slice_layers = 0;

			if (max_clipping_plane - min_clipping_plane != 0)
				all_calculated_value[i].value_of_less_clipping_plane = 1 - double(Tree_nodes_cut_layers[candidate_nodes[i]].size() - min_clipping_plane) / double(max_clipping_plane - min_clipping_plane);
			else
				all_calculated_value[i].value_of_less_clipping_plane = 0;
		}*/
		/////////////////


		//Calculate the accumulated value//
	int* sum_area_S = new int[all_calculated_value.size()];
	double* sum_slice_layers = new double[all_calculated_value.size()];
	int* sum_clipping_plane = new int[all_calculated_value.size()];
	double* sum_larger_base = new double[all_calculated_value.size()];
	double* sum_covering_points = new double[all_calculated_value.size()];
	double* sum_fragile = new double[all_calculated_value.size()];
	double* sum_orientation = new double[all_calculated_value.size()];
	double* sum_projected = new double[all_calculated_value.size()];
	for (int i = 0; i < all_calculated_value.size(); i++) {
		sum_area_S[i] = 0;
		sum_slice_layers[i] = 0;
		sum_clipping_plane[i] = 0;
		sum_larger_base[i] = 0;
		sum_covering_points[i] = 0;
		sum_orientation[i] = 0;
		sum_fragile[i] = all_calculated_value[i].value_of_fragile;
		//cout << "((((((" << sum_fragile[i] << endl;
		sum_projected[i] = all_calculated_value[i].value_of_projected;
		int index_node = candidate_nodes[i];
		pure_value[i].value_of_orientation = save_ori[index_node][0].dot(save_ori[pre_tree_nodes[index_node]][0]);
		pure_value[i].value_of_fragile = all_calculated_value[i].value_of_fragile;
		pure_value[i].value_of_projected = all_calculated_value[i].value_of_projected;
		while (pre_tree_nodes[index_node] != -1) {
			sum_area_S[i] += int(final_pathes_include_S[index_node].size());
			double temp_sum = 0;
			for (int j = 0; j < Tree_nodes[index_node].size(); j++) {
				for (int k = 0; k < Tree_nodes[index_node][j].size() - 1; k++) {
					sum_slice_layers[i] += distance3d(Tree_nodes[index_node][j][k], Tree_nodes[index_node][j][k + 1]);
				}
			}
			sum_clipping_plane[i] += int(Tree_nodes_cut_layers[index_node].size());
			sum_larger_base[i] += Tree_nodes_larger_base[index_node];

			//sum_covering_points[i] += int(final_pathes_include_covering_points[index_node].size());
			for (int j = 0; j < final_pathes_include_covering_points[index_node].size(); j++)
				sum_covering_points[i] += int(all_the_covering_points[final_pathes_include_covering_points[index_node][j]].size());
			index_node = pre_tree_nodes[index_node];
		}
		index_node = candidate_nodes[i];
		while (pre_tree_nodes[index_node] != 0) {
			sum_orientation[i] += save_ori[index_node][0].dot(save_ori[pre_tree_nodes[index_node]][0]);
			index_node = pre_tree_nodes[index_node];
		}
	}
	////////////////////////////////////


	//normalization - accumulated//
	int max_area_S = -MAX_I, min_area_S = MAX_I;
	double max_slice_layers = -MAX_I, min_slice_layers = MAX_I;
	int max_clipping_plane = -MAX_I, min_clipping_plane = MAX_I;
	double max_large_base = -MAX_I, min_large_base = MAX_I;
	double max_covering_points = -MAX_I, min_covering_points = MAX_I;
	double max_fragile = -MAX_I, min_fragile = MAX_I;
	double max_orientation = -MAX_I, min_orientation = MAX_I;
	double max_projected = -MAX_I, min_projected = MAX_I;
	for (int i = 0; i < all_calculated_value.size(); i++) {
		max_large_base = max(max_large_base, sum_larger_base[i]);
		min_large_base = min(min_large_base, sum_larger_base[i]);

		max_area_S = max(max_area_S, sum_area_S[i]);
		min_area_S = min(min_area_S, sum_area_S[i]);

		max_slice_layers = max(max_slice_layers, sum_slice_layers[i]);
		min_slice_layers = min(min_slice_layers, sum_slice_layers[i]);

		max_clipping_plane = max(max_clipping_plane, sum_clipping_plane[i]);
		min_clipping_plane = min(min_clipping_plane, sum_clipping_plane[i]);

		max_covering_points = max(max_covering_points, sum_covering_points[i]);
		min_covering_points = min(min_covering_points, sum_covering_points[i]);

		max_fragile = max(max_fragile, sum_fragile[i]);
		min_fragile = min(min_fragile, sum_fragile[i]);

		max_orientation = max(max_orientation, sum_orientation[i]);
		min_orientation = min(min_orientation, sum_orientation[i]);

		max_projected = max(max_projected, sum_projected[i]);
		min_projected = min(min_projected, sum_projected[i]);
	}

	for (int i = 0; i < all_calculated_value.size(); i++) {
		if (max_large_base - min_large_base != 0)
			all_calculated_value[i].large_base = (sum_larger_base[i] - min_large_base) / (max_large_base - min_large_base);
		else
			all_calculated_value[i].large_base = 1;  //全部暂时改为1，之后还是需要-1用于标识

		if (max_area_S - min_area_S != 0)
			all_calculated_value[i].value_of_area_S = 1 - double(sum_area_S[i] - min_area_S) / double(max_area_S - min_area_S);
		else
			all_calculated_value[i].value_of_area_S = 1;

		if (max_slice_layers - min_slice_layers != 0)
			all_calculated_value[i].value_of_more_slice_layers = double(sum_slice_layers[i] - min_slice_layers) / double(max_slice_layers - min_slice_layers);
		else
			all_calculated_value[i].value_of_more_slice_layers = 1;

		if (max_clipping_plane - min_clipping_plane != 0)
			all_calculated_value[i].value_of_less_clipping_plane = 1 - double(sum_clipping_plane[i] - min_clipping_plane) / double(max_clipping_plane - min_clipping_plane);
		else
			all_calculated_value[i].value_of_less_clipping_plane = 1;

		if (max_covering_points - min_covering_points != 0)
			all_calculated_value[i].value_of_covering_points = double(sum_covering_points[i] - min_covering_points) / double(max_covering_points - min_covering_points);
		else
			all_calculated_value[i].value_of_covering_points = 1;

		if (max_fragile - min_fragile != 0)
			all_calculated_value[i].value_of_fragile = 1 - double(sum_fragile[i] - min_fragile) / double(max_fragile - min_fragile);
		else
			all_calculated_value[i].value_of_fragile = 1;

		if (max_orientation - min_orientation != 0)
			all_calculated_value[i].value_of_orientation = double(sum_orientation[i] - min_orientation) / double(max_orientation - min_orientation);
		else
			all_calculated_value[i].value_of_orientation = 1;

		/*if (height_of_beam_search == 3) {
			if (max_projected - min_projected != 0)
				all_calculated_value[i].value_of_projected = double(sum_projected[i] - min_projected) / double(max_projected - min_projected);
			else
				all_calculated_value[i].value_of_projected = -1;
		}*/
		//else {
		if (max_projected - min_projected != 0)
			all_calculated_value[i].value_of_projected = 1 - double(sum_projected[i] - min_projected) / double(max_projected - min_projected);
		else
			all_calculated_value[i].value_of_projected = 1;
		//}
	}
	/////////////////

	//sort
	for (int i = 0; i < candidate_nodes.size(); i++) {
		for (int j = i + 1; j < candidate_nodes.size(); j++) {
			if ((all_calculated_value[i].value_of_area_S * W_less_area_S + all_calculated_value[i].value_of_more_slice_layers * W_more_slice_layers + all_calculated_value[i].value_of_less_clipping_plane * W_less_clipping_plane + all_calculated_value[i].large_base * W_larger_base + all_calculated_value[i].value_of_covering_points * W_covering_points + all_calculated_value[i].value_of_fragile * W_fragile + all_calculated_value[i].value_of_orientation * W_orientation + all_calculated_value[i].value_of_projected * W_projected)
				< (all_calculated_value[j].value_of_area_S * W_less_area_S + all_calculated_value[j].value_of_more_slice_layers * W_more_slice_layers + all_calculated_value[j].value_of_less_clipping_plane * W_less_clipping_plane + all_calculated_value[j].large_base * W_larger_base + all_calculated_value[j].value_of_covering_points * W_covering_points + all_calculated_value[j].value_of_fragile * W_fragile + all_calculated_value[j].value_of_orientation * W_orientation + all_calculated_value[j].value_of_projected * W_projected))
			{
				swap(candidate_nodes[i], candidate_nodes[j]);
				swap(all_calculated_value[i], all_calculated_value[j]);
				swap(pure_value[i], pure_value[j]);
			}
		}
	}


	//if (height_of_beam_search == 6) {
		//cout << endl << "&&" << all_calculated_value[0].number_of_remaining_face << endl;
		//cout << endl << "&&" << all_calculated_value[1].number_of_remaining_face << endl;
		//}

	//先删除
	cout << "#####################" << height_of_beam_search << " " << id_continue << "#####################" << endl;
	//if (height_of_beam_search == 2) {  //julia_vase
	//	int cont = 0;
	//	while (cont < 2) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	// 
	//if (height_of_beam_search == 3) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//else if (height_of_beam_search == 5) {
	//	int cont = 0;
	//	while (cont <1) {  //6
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 2 && id_continue == 0) {
	//	int cont = 0;
	//	while (cont < 1) {  //6
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}

	//if (height_of_beam_search == 4) {
	//	int cont = 0;
	//	while (cont < 5) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 5) {
	//	int cont = 0;
	//	while (cont < 8) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 6) {
	//	int cont = 0;
	//	while (cont < 13) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 7) {
	//	int cont = 0;
	//	while (cont < 0) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 13) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 14) {
	//	int cont = 0;
	//	while (cont < 2) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
	//if (height_of_beam_search == 18) {
	//	int cont = 0;
	//	while (cont < 1) {  //5
	//		candidate_nodes.erase(candidate_nodes.begin());
	//		all_calculated_value.erase(all_calculated_value.begin());
	//		cont++;
	//	}
	//}
}

void HybridManufacturing::sort_candidate_nodes(vector<int>& candidate_nodes, vector<vector<int>> Tree_nodes_for_S)
{
	/*for (int i = 0; i < candidate_nodes.size(); i++) {
		for (int j = i + 1; j < candidate_nodes.size(); j++) {
			if (Tree_nodes_for_S[candidate_nodes[i]].size() < Tree_nodes_for_S[candidate_nodes[j]].size()) {
				swap(candidate_nodes[i], candidate_nodes[j]);
			}
		}
	}*/
}

void HybridManufacturing::subtractive_remove_output(const vector<TRiangle>& need_detect_triangle, const Slicer_2& current_slicer, int height_of_beam_search)
{
	std::string filename = ".\\vis\\block_patch-" + to_string(height_of_beam_search) + "_" + ".obj";
	std::ofstream file(filename);
	if (!file)
	{
		std::cout << "subtractive_remove_output !file" << std::endl;
	}

	for (auto& p : current_slicer.positions) {
		file << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
	}
	for (auto t : need_detect_triangle) {
		for (int i = 0; i < 3; ++i) t[i]++;
		file << "f " << t[0] << " " << t[1] << " " << t[2] << "\n";
	}
}

bool HybridManufacturing::CheckToolCollisionWithCell(
	const Eigen::Vector3d& center_point,
	const std::vector<Eigen::MatrixXd>& target_cell_vertices,
	double max_z_target,
	const cutter& nozzle,
	double z_threshold_divisor,
	double xy_tolerance) const
{
	(void)xy_tolerance;

	const double z_threshold = nozzle.cylinder_r / z_threshold_divisor;
	const double height_diff = max_z_target - center_point.z();

	// 快速排除：目标在刀尖球下方
	if (height_diff <= z_threshold) {
		return false;
	}

	// 快速确认：目标超出刀具最大高度
	if (height_diff > nozzle.total_height) {
		return true;
	}

	// 使用第一个顶点进行粗筛
	if (!target_cell_vertices.empty()) {
		const double dx = target_cell_vertices[0](0, 0) - center_point.x();
		const double dy = target_cell_vertices[0](1, 0) - center_point.y();
		const double dist_xy_sq = dx * dx + dy * dy;

		if (height_diff > nozzle.cylinder_height_threshold) {
			if (dist_xy_sq > nozzle.carriage_check_radius_sq) {
				return false;
			}
		}
		else {
			if (dist_xy_sq > nozzle.cylinder_check_radius_sq) {
				return false;
			}
		}
	}

	// 精确检测每个边界顶点
	for (const auto& vertex : target_cell_vertices) {
		const double vx = vertex(0, 0);
		const double vy = vertex(1, 0);
		const double vz = vertex(2, 0);
		const double diff_z = vz - center_point.z();

		if (diff_z <= z_threshold) {
			continue;
		}

		const double dx = vx - center_point.x();
		const double dy = vy - center_point.y();
		const double dist_xy_sq = dx * dx + dy * dy;

		if (diff_z <= nozzle.cylinder_height_threshold) { //nozzle.cylinder_r + nozzle.cylinder_height
			if (dist_xy_sq < nozzle.cylinder_r_sq) {
				return true;
			}
		}
		else if (diff_z <= nozzle.total_height) { //nozzle.cylinder_r + nozzle.cylinder_height + nozzle.carriage_height
			if (dist_xy_sq < nozzle.carriage_r_sq) {
				return true;
			}
		}
	}

	return false;
}
