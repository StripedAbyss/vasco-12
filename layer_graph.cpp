#include "layer_graph.h"

namespace {
static bool SurfaceMeshContourIsSelfSupport(
	const std::vector<int>& face_ids,
	const std::vector<Eigen::Vector3d>& face_normals)
{
	if (face_ids.empty()) {
		return true;
	}
	const Eigen::Vector3d base_normal(0, 0, 1);
	int cont_not_set_support = 0;
	for (int face_id : face_ids) {
		if (face_id < 0 || face_id >= static_cast<int>(face_normals.size())) {
			continue;
		}
		const Eigen::Vector3d& face_normal = face_normals[face_id];
		if (face_normal.dot(base_normal) + sin(PI / 3.6) < 0) {
			++cont_not_set_support;
		}
	}
	return double(cont_not_set_support) / double(face_ids.size()) <= 0.2;
}

static std::vector<Eigen::Vector2d> SliceContourToPolygon(const Polyline_type& contour)
{
	std::vector<Eigen::Vector2d> polygon;
	polygon.reserve(contour.size());
	for (const auto& p : contour) {
		polygon.emplace_back(p.x(), p.y());
	}
	return polygon;
}

// Tests whether a point lies in the offset band of any component boundary.
// Both the outer contour and every hole contour contribute supporting walls.
static bool PointInsideComponentBoundaryBand(
	Polygon& outer_boundary_band,
	std::vector<Polygon>& hole_boundary_bands,
	const Eigen::Vector2d& point)
{
	if (outer_boundary_band.JudgePointInside(point)) {
		return true;
	}
	for (Polygon& hole_boundary_band : hole_boundary_bands) {
		if (hole_boundary_band.JudgePointInside(point)) {
			return true;
		}
	}
	return false;
}

// Tests whether any outer or hole boundary pair is closer than the nozzle radius.
// An AABB rejection keeps distant component pairs out of the quadratic point test.
static bool ComponentBoundariesWithinDistance(
	const Data& data,
	std::size_t first_layer,
	std::size_t first_component,
	std::size_t second_layer,
	std::size_t second_component,
	double distance_squared)
{
	const auto test_boundary_pair = [distance_squared](
		const std::vector<Eigen::Vector2d>& first,
		const std::vector<Eigen::Vector2d>& second) {
		if (first.empty() || second.empty()) {
			return false;
		}

		double first_min_x = first.front().x();
		double first_max_x = first.front().x();
		double first_min_y = first.front().y();
		double first_max_y = first.front().y();
		for (const Eigen::Vector2d& point : first) {
			first_min_x = std::min(first_min_x, point.x());
			first_max_x = std::max(first_max_x, point.x());
			first_min_y = std::min(first_min_y, point.y());
			first_max_y = std::max(first_max_y, point.y());
		}

		double second_min_x = second.front().x();
		double second_max_x = second.front().x();
		double second_min_y = second.front().y();
		double second_max_y = second.front().y();
		for (const Eigen::Vector2d& point : second) {
			second_min_x = std::min(second_min_x, point.x());
			second_max_x = std::max(second_max_x, point.x());
			second_min_y = std::min(second_min_y, point.y());
			second_max_y = std::max(second_max_y, point.y());
		}

		const double distance = std::sqrt(std::max(0.0, distance_squared));
		if (first_max_x + distance < second_min_x
			|| second_max_x + distance < first_min_x
			|| first_max_y + distance < second_min_y
			|| second_max_y + distance < first_min_y) {
			return false;
		}

		for (const Eigen::Vector2d& a : first) {
			for (const Eigen::Vector2d& b : second) {
				if ((a - b).squaredNorm() < distance_squared) {
					return true;
				}
			}
		}
		return false;
		};

	std::vector<const std::vector<Eigen::Vector2d>*> first_boundaries{
		&data.slice_points[first_layer][first_component]
	};
	std::vector<const std::vector<Eigen::Vector2d>*> second_boundaries{
		&data.slice_points[second_layer][second_component]
	};
	for (const auto& hole : data.slice_points_holes[first_layer][first_component]) {
		first_boundaries.push_back(&hole);
	}
	for (const auto& hole : data.slice_points_holes[second_layer][second_component]) {
		second_boundaries.push_back(&hole);
	}

	for (const auto* first : first_boundaries) {
		for (const auto* second : second_boundaries) {
			if (test_boundary_pair(*first, *second)) {
				return true;
			}
		}
	}
	return false;
}
}

Layer_Graph::Layer_Graph(const Data& data)
{
	this->total_node_num = data.total_node_num;
	this->in_degree.resize(this->total_node_num, 0);
	this->out_degree.resize(this->total_node_num, 0);
	this->node_visited.resize(this->total_node_num, false);
	this->edge_deleted.resize(this->total_node_num * this->total_node_num, false);
	this->data = data;
	this->G.resize(maxn);
	this->G_2.resize(maxn);
	this->G_3.resize(maxn);
	this->cont_normal_dependency_edges = 0;
}

Layer_Graph::~Layer_Graph()
{
}


void Layer_Graph::creat_ball(string file_name, std::vector<Eigen::MatrixXd> vis_points)
{
	Eigen::MatrixXd V_2;
	Eigen::MatrixXi F_2;

	igl::readOBJ("ball.obj", V_2, F_2);
	ofstream all_balls(".\\vis\\" + file_name + "_unaccessivle_points.obj");
	for (int i = 0; i < vis_points.size(); i++) {
		for (int j = 0; j < V_2.rows(); j++)
			all_balls << "v " << V_2(j, 0) + vis_points[i](0,0) << " " << V_2(j, 1) + vis_points[i](1, 0) << " " << V_2(j, 2) + vis_points[i](2, 0)<< " 0.9"<< " 0.05" <<" 0.05" << endl;
		for (int j = 0; j < F_2.rows(); j++)
			all_balls << "f " << F_2(j, 0) + i * V_2.rows() + 1 << " " << F_2(j, 1) + i * V_2.rows() + 1 << " " << F_2(j, 2) + i * V_2.rows() + 1 << endl;
	}
}

void Layer_Graph::GetTrianglesForSurfaceMeshSlices(
	const std::vector<SurfaceMeshSliceData>& slices,
	const Eigen::Vector3d& vectorAfter,
	int height_of_beam_search,
	int id_continue)
{
	std::vector<std::vector<std::vector<int>>> contour_face_ids;
	std::vector<Eigen::Vector3d> face_normals;
	for (std::size_t layer_id = 0; layer_id < data.slice_points.size(); ++layer_id) {
		contour_face_ids.push_back({});
		if (layer_id >= slices.size()) {
			continue;
		}
		const auto& slice = slices[layer_id];
		for (std::size_t component_id = 0;
			component_id < data.slice_points[layer_id].size();
			++component_id) {
			contour_face_ids.back().push_back({});
			auto& component_face_ids = contour_face_ids.back().back();
			if (layer_id >= data.source_contour_ids.size()
				|| layer_id >= data.source_hole_contour_ids.size()
				|| component_id >= data.source_contour_ids[layer_id].size()
				|| component_id >= data.source_hole_contour_ids[layer_id].size()) {
				continue;
			}
			const int source_id = data.source_contour_ids[layer_id][component_id];
			if (source_id >= 0
				&& static_cast<std::size_t>(source_id) < slice.contour_face_ids.size()) {
				for (const auto& face_index : slice.contour_face_ids[static_cast<std::size_t>(source_id)]) {
					component_face_ids.push_back(static_cast<int>(face_index));
				}
			}
			for (int hole_source_id : data.source_hole_contour_ids[layer_id][component_id]) {
				if (hole_source_id < 0
					|| static_cast<std::size_t>(hole_source_id) >= slice.contour_face_ids.size()) {
					continue;
				}
				for (const auto& face_index : slice.contour_face_ids[static_cast<std::size_t>(hole_source_id)]) {
					component_face_ids.push_back(static_cast<int>(face_index));
				}
			}
		}
		for (const auto& normal : slice.face_normals) {
			face_normals.push_back(normal);
		}
	}
	GetTrianglesForLayersFromMesh(contour_face_ids, face_normals, vectorAfter, height_of_beam_search, id_continue);
}

void Layer_Graph::GenerateDependencyEdgesFromSurfaceMeshSlices( //������ֻ���ݵ��Ƿ�����һ����������������������ߣ����Ըĳɸ��ݱ��Ƿ񴩹���һ���������������������
	const std::vector<SurfaceMeshSliceData>& slices)
{
	for (std::size_t i = 1; i < data.slice_points.size(); ++i) {
		std::vector<Polygon> last_layer_polygons;
		std::vector<std::vector<Polygon>> last_layer_hole_polygons;
		last_layer_polygons.reserve(data.slice_points[i - 1].size());
		last_layer_hole_polygons.resize(data.slice_points[i - 1].size());
		for (std::size_t component_id = 0;
			component_id < data.slice_points[i - 1].size();
			++component_id) {
			last_layer_polygons.emplace_back(
				ConstructPolygonPoints(data.slice_points[i - 1][component_id], dependence_offset));
			for (const auto& hole : data.slice_points_holes[i - 1][component_id]) {
				last_layer_hole_polygons[component_id].emplace_back(
					ConstructPolygonPoints(hole, dependence_offset));
			}
		}

		for (std::size_t j = 0; j < data.slice_points[i].size(); ++j) {
			for (std::size_t m = 0; m < last_layer_polygons.size(); ++m) {
				for (const Eigen::Vector2d& test_point : data.slice_points[i][j]) {
					if (!PointInsideComponentBoundaryBand(
						last_layer_polygons[m],
						last_layer_hole_polygons[m],
						test_point)) {
						continue;
					}
					const int from_id = data.index_inv[std::make_pair(
						static_cast<int>(i - 1),
						static_cast<int>(m))];
					const int to_id = data.index_inv[std::make_pair(
						static_cast<int>(i),
						static_cast<int>(j))];
					this->AddEdge(from_id, to_id);
					temp_edges.push_back(make_pair(from_id, to_id));
					break;
				}
			}
		}
	}
}

void Layer_Graph::CollisionDetectionForAdditiveManufacturingFromSurfaceMeshSlices(
	const std::vector<SurfaceMeshSliceData>& slices,
	nozzle the_nozzle)
{
	if (slices.empty()) {
		std::cout << "[Layer_Graph::CollisionDetectionForAdditiveManufacturingFromSurfaceMeshSlices] Warning: No slices provided. Collision detection skipped." << std::endl;
		return;
	}


	vector<pair<int, int>> temp_collision_edges;

	for (std::size_t i = 0; i < data.slice_points.size(); ++i) {
		for (std::size_t j = 0; j < data.slice_points[i].size(); ++j) {
			for (std::size_t ii = i + 1; ii < data.slice_points.size(); ++ii) {
				double circle_r;
				if (slices[ii].layer_z - slices[i].layer_z < 0) {
					std::cout << "[Layer_Graph::CollisionDetectionForAdditiveManufacturingFromSurfaceMeshSlices] Warning: Layer " << ii << " is below layer " << i << ". Skipping collision detection for this pair of layers." << std::endl;
					continue;
				}
				double layer_distance = slices[ii].layer_z - slices[i].layer_z;
				if (layer_distance < the_nozzle.nozzle_H_half) {
					circle_r = the_nozzle.lowwer_surface_r + layer_distance * (the_nozzle.upper_surface_r - the_nozzle.lowwer_surface_r) / the_nozzle.nozzle_H_half;
				}
				//if ((ii - i) * dh < the_nozzle.nozzle_H_half)
				//	circle_r = the_nozzle.lowwer_surface_r + (ii - i) * dh * (the_nozzle.upper_surface_r - the_nozzle.lowwer_surface_r) / the_nozzle.nozzle_H_half;
				else {
					circle_r = the_nozzle.upper_surface_r;
				}
				for (std::size_t jj = 0; jj < data.slice_points[ii].size(); ++jj) {
					const int lower_layer_id = static_cast<int>(i);
					const int lower_component_id = static_cast<int>(j);
					const int upper_layer_id = static_cast<int>(ii);
					const int upper_component_id = static_cast<int>(jj);
					bool jud_collision = false;
					if (layer_distance > the_nozzle.nozzle__H_total) {
						jud_collision = true;
						temp_collision_edges.push_back(make_pair(
							data.index_inv[std::make_pair(lower_layer_id, lower_component_id)],
							data.index_inv[std::make_pair(upper_layer_id, upper_component_id)]));
						continue;
					}
					else {
						jud_collision = ComponentBoundariesWithinDistance(
							data,
							i,
							j,
							ii,
							jj,
							circle_r * circle_r);
					}
					if (jud_collision) {
						const int from_id =
							data.index_inv[std::make_pair(lower_layer_id, lower_component_id)];
						const int to_id =
							data.index_inv[std::make_pair(upper_layer_id, upper_component_id)];
						//this->AddEdge_2(from_id, to_id);
						//temp_edges.push_back(make_pair(from_id, to_id));
						temp_collision_edges.push_back(make_pair(from_id, to_id));
					}
				}
			}
		}
	}


	// Build the dependency matrix without raw pointers.
	std::vector<std::vector<bool>> dependency_relationship(
		data.total_node_num,
		std::vector<bool>(data.total_node_num, false));

	for (int layer_id = 1; layer_id < static_cast<int>(data.slice_points.size()); ++layer_id) {
		for (int point_id = 0; point_id < static_cast<int>(data.slice_points[layer_id].size()); ++point_id) {
			const int target_node = data.index_inv[std::make_pair(layer_id, point_id)];
			std::vector<int> current_frontier{ target_node };

			for (int remaining_layers = layer_id; remaining_layers > 0 && !current_frontier.empty(); --remaining_layers) {
				std::vector<int> next_frontier;
				for (const auto& edge : temp_edges) {
					if (std::find(current_frontier.begin(), current_frontier.end(), edge.second) == current_frontier.end()) {
						continue;
					}

					dependency_relationship[edge.first][target_node] = true;
					next_frontier.push_back(edge.first);
				}
				current_frontier = std::move(next_frontier);
			}
		}
	}

	cont_normal_dependency_edges = temp_edges.size();


	//cout << "the number of dependency edges: " << temp_edges.size() << endl;
	for (const auto& collision_edge : temp_collision_edges) {
		const int from = collision_edge.first;
		const int to = collision_edge.second;
		if (dependency_relationship[from][to]) {
			continue;
		}

		dependency_relationship[from][to] = true;
		for (int node_id = 0; node_id < data.total_node_num; ++node_id) {
			if (dependency_relationship[to][node_id]) {
				dependency_relationship[from][node_id] = true;
			}
		}

		this->AddEdge_2(from, to);
		temp_edges.push_back(make_pair(from, to));
	}


	//cout << endl << "time..." << double(end_time_9 - start_time_9) / CLOCKS_PER_SEC << endl;
}

void Layer_Graph::BuildLayerGraphFromSurfaceMeshSlices(
	const std::vector<SurfaceMeshSliceData>& slices,
	nozzle the_nozzle)
{
	this->total_node_num = data.total_node_num;
	this->in_degree.assign(this->total_node_num, 0);
	this->out_degree.assign(this->total_node_num, 0);
	this->node_visited.assign(this->total_node_num, false);
	this->edge_deleted.assign(this->total_node_num * this->total_node_num, false);
	this->G.resize(maxn);
	this->G_2.resize(maxn);
	this->G_3.resize(maxn);
	this->temp_edges.clear();
	this->cont_normal_dependency_edges = 0;
	this->all_triangles_of_layers.clear();
	this->is_the_layer_self_suppot.clear();
	this->GetTrianglesForSurfaceMeshSlices(slices, Eigen::Vector3d(0, 0, 1), 0, 0);
	this->GenerateDependencyEdgesFromSurfaceMeshSlices(slices);
	this->CollisionDetectionForAdditiveManufacturingFromSurfaceMeshSlices(slices, the_nozzle);
}

void Layer_Graph::GetTrianglesForLayersFromMesh(
	const std::vector<std::vector<std::vector<int>>>& contour_face_ids,
	const std::vector<Eigen::Vector3d>& face_normals,
	const Eigen::Vector3d& vectorAfter,
	int height_of_beam_search,
	int id_continue)
{
	is_the_layer_self_suppot.resize(total_node_num);
	int cont_num = 0;
	all_triangles_of_layers.resize(total_node_num);
	const Eigen::Vector3d base_normal(0, 0, 1);

	for (int i = 0; i < contour_face_ids.size(); i++) {
		for (int j = 0; j < contour_face_ids[i].size(); j++) {
			is_the_layer_self_suppot[cont_num] = true;
			int cont_not_set_support = 0;
			const auto& face_ids = contour_face_ids[i][j];
			for (int k = 0; k < face_ids.size(); k++) {
				int face_id = face_ids[k];
				if (face_id < 0 || face_id >= face_normals.size()) {
					continue;
				}
				const Eigen::Vector3d& face_normal = face_normals[face_id];
				bool jud_self_support = (face_normal.dot(base_normal) + sin(PI / 3.6) >= 0);
				//bool jud_self_support = (face_normal.dot(base_normal) + sin(PI / 4) >= 0);
				if (!jud_self_support) {
					cont_not_set_support++;
				}
			}
			if (!face_ids.empty()) {
				if (double(cont_not_set_support) / double(face_ids.size()) > 0.2) {
					is_the_layer_self_suppot[cont_num] = false;
				}
			}
			cont_num++;
		}
	}
}


void Layer_Graph::GetTrianglesForLayers(vector<vector<vector<Vertex>>> all_slice_points, std::vector<map<pair<Vertex, Vertex>, Triangle*>> map_segment_triangles, vector<Vertex> all_vertex, Eigen::Vector3d vectorAfter, int height_of_beam_search, int id_continue)
{	is_the_layer_self_suppot.resize(total_node_num);
	int cont_num = 0;
	all_triangles_of_layers.resize(total_node_num);
	//cout << "%%" << data.slice_points.size() << endl;
	for (int i = 0; i < data.slice_points.size(); i++) {
		for (int j = 0; j < data.slice_points[i].size(); j++) {
			is_the_layer_self_suppot[cont_num] = true;
			int cont_not_set_support = 0;
			for (int k = 1; k < data.slice_points[i][j].size(); k++) {
				Vertex v_left = all_slice_points[i][j][k - 1];
				Vertex v_right = all_slice_points[i][j][k];
				pair<Vertex, Vertex> temp_pair(v_left, v_right);
				if (map_segment_triangles[i].count(temp_pair) <= 0)
					continue;
				Triangle* tri = map_segment_triangles[i][temp_pair];
				if (tri == nullptr) {
					std::cout << "[Layer_Graph::GetTrianglesForLayers] Triangle is null" << std::endl;
					continue;
				}
				all_triangles_of_layers[cont_num].push_back(tri);
				//cout << "d" << endl;
				Vertex* v1 = tri->vertices_2[0];
				//cout << "k" << endl;
				Vertex* v2 = tri->vertices_2[1];
				//cout << "m" << endl;
				Vertex* v3 = tri->vertices_2[2];
				//cout << "d" << endl;
				//double ans = (v2->x - v1->x) * (v2->y - v3->y) - (v2->y - v1->y) * (v2->x - v3->x);
				//if (ans > 0)	//is clockwise
				//	swap(v2, v3);
				double na = (v2->y - v1->y) * (v3->z - v1->z) - (v2->z - v1->z) * (v3->y - v1->y);
				double nb = (v2->z - v1->z) * (v3->x - v1->x) - (v2->x - v1->x) * (v3->z - v1->z);
				double nc = (v2->x - v1->x) * (v3->y - v1->y) - (v2->y - v1->y) * (v3->x - v1->x);
				Eigen::Vector3d face_normal(na, nb, nc);
				face_normal.normalize();
				Eigen::Vector3d base_normal (0, 0, 1);
				base_normal.normalize();
				bool jud_self_support;
				//cout << "#####" << height_of_beam_search << " " << id_continue << endl;
				//if (height_of_beam_search != 6 || id_continue != 0)
				//	jud_self_support = (face_normal.dot(base_normal) + sin(PI / 10) >= 0);   //PI / 3.6
				//else
					jud_self_support = (face_normal.dot(base_normal) + sin(PI / 3.6) >= 0);   //PI / 3.6
				if (jud_self_support == false) {
					//cont_not_set_support++;   //ע��ʱ�ر�
				}
				//cout << "%%" << height_of_beam_search << " " << id_continue << endl;
				//if (height_of_beam_search <= 2) {
					if (double(cont_not_set_support) / double(data.slice_points[i][j].size()) > 0.2) {    //data.slice_points[i][j].size()) > 0.2 
						is_the_layer_self_suppot[cont_num] = false;
						break;
					}
				//}
				//else {
				//	if (double(cont_not_set_support) / double(data.slice_points[i][j].size()) > 1.0) {    //data.slice_points[i][j].size()) > 0.2  
				//		is_the_layer_self_suppot[cont_num] = false;
				//		break;
				//	}
				//}
			}
			cont_num++;
		}
	}
}


void Layer_Graph::GenerateDependencyEdges()
{
	for (std::size_t layer_id = 1; layer_id < data.slice_points.size(); ++layer_id) {
		std::vector<Polygon> lower_outer_polygons;
		std::vector<std::vector<Polygon>> lower_hole_polygons;
		lower_outer_polygons.reserve(data.slice_points[layer_id - 1].size());
		lower_hole_polygons.resize(data.slice_points[layer_id - 1].size());

		for (std::size_t component_id = 0;
			component_id < data.slice_points[layer_id - 1].size();
			++component_id) {
			lower_outer_polygons.emplace_back(
				ConstructPolygonPoints(
					data.slice_points[layer_id - 1][component_id],
					dependence_offset));
			for (const auto& hole : data.slice_points_holes[layer_id - 1][component_id]) {
				lower_hole_polygons[component_id].emplace_back(
					ConstructPolygonPoints(hole, dependence_offset));
			}
		}

		for (std::size_t upper_component_id = 0;
			upper_component_id < data.slice_points[layer_id].size();
			++upper_component_id) {
			for (std::size_t lower_component_id = 0;
				lower_component_id < lower_outer_polygons.size();
				++lower_component_id) {
				for (const Eigen::Vector2d& test_point : data.slice_points[layer_id][upper_component_id]) {
					if (!PointInsideComponentBoundaryBand(
						lower_outer_polygons[lower_component_id],
						lower_hole_polygons[lower_component_id],
						test_point)) {
						continue;
					}

					const int from_id = data.index_inv[std::make_pair(
						static_cast<int>(layer_id - 1),
						static_cast<int>(lower_component_id))];
					const int to_id = data.index_inv[std::make_pair(
						static_cast<int>(layer_id),
						static_cast<int>(upper_component_id))];
					this->AddEdge(from_id, to_id);
					temp_edges.emplace_back(from_id, to_id);
					break;
				}
			}
		}
	}
}

void Layer_Graph::BuildLayerGraph(nozzle the_nozzle)
{
	clock_t start_time, end_time;
	start_time = clock();
    vector<pair<int, int>> temp_collision_edges;
	vector<pair<int, int>> all_id_layers;
	Eigen::Vector2d dir;
	for (int i = 1; i < data.slice_points.size(); i++) {
		//get (i-1) layer segment's contour
		std::vector<Polygon> Last_layer_polygons;
		for (int j = 0; j < data.slice_points[i - 1].size(); j++) {
			Last_layer_polygons.push_back(Polygon(ConstructPolygonPoints(data.slice_points[i - 1][j], dependence_offset)));
		}
		// build graph
		for (int j = 0; j < data.slice_points[i].size(); j++) {
			for (int m = 0; m < Last_layer_polygons.size(); m++) {
				for (int k = 0; k < data.slice_points[i][j].size(); k++) {
					if (Last_layer_polygons[m].JudgePointInside(data.slice_points[i][j][k])) {
						this->AddEdge(data.index_inv[std::make_pair(i - 1, m)], data.index_inv[std::make_pair(i, j)]);
						temp_edges.push_back(make_pair(data.index_inv[std::make_pair(i - 1, m)], data.index_inv[std::make_pair(i, j)]));
						all_id_layers.push_back(make_pair(i - 1 , i));
						break;
					}
				}
			}
		}
	}
	cont_normal_dependency_edges = temp_edges.size();

	//collision detect, add dependency edge
	for (int i = 0; i < data.slice_points.size(); i++) {
		for (int j = 0; j < data.slice_points[i].size(); j++) {
			for (int ii = i + 1; ii < data.slice_points.size(); ii++) {
				double circle_r;
				if ((ii - i) * dh < the_nozzle.nozzle_H_half)
					circle_r = the_nozzle.lowwer_surface_r + (ii - i) * dh * (the_nozzle.upper_surface_r - the_nozzle.lowwer_surface_r) / the_nozzle.nozzle_H_half;
				else
					circle_r = the_nozzle.upper_surface_r;

				for (int jj = 0; jj < data.slice_points[ii].size(); jj++) {
					bool jud_collision = false;
					if ((ii - i) * dh > the_nozzle.nozzle__H_total) //exceed nozzle_H
					{
						jud_collision = true;
						temp_collision_edges.push_back(make_pair(data.index_inv[std::make_pair(i, j)], data.index_inv[std::make_pair(ii, jj)]));
						continue;
					}
					for (int k = 0; k < data.slice_points[i][j].size(); k+=20) { //step == 2
						for (int kk = 0; kk < data.slice_points[ii][jj].size(); kk+=20) {
                      if (pow(data.slice_points[ii][jj][kk].x() - data.slice_points[i][j][k].x(),2) + pow(data.slice_points[ii][jj][kk].y() - data.slice_points[i][j][k].y(),2) - pow(circle_r,2) < 0) {
								jud_collision = true;
								break;
							}
						}
						if (jud_collision == true)
							break;
					}
					if (jud_collision == true) {
						temp_collision_edges.push_back(make_pair(data.index_inv[std::make_pair(i, j)], data.index_inv[std::make_pair(ii, jj)]));
					}
				}
			}
		}
	}
	
	//establish hase graph
	bool** dependency_relationship = new bool*[10000];
	for (int i = 0;i < 10000;i++)
		dependency_relationship[i] = new bool[10000];
	for (int i = 0;i < 10000;i++)
		for (int j = 0;j < 10000;j++)
			dependency_relationship[i][j] = false;

	for (int i = 1; i < data.slice_points.size(); i++) {
		for (int j = 0; j < data.slice_points[i].size(); j++) {
			int id_layer = i;
			vector<int> id_current_layer;
			vector<int> temp_id_current_layer;
			id_current_layer.push_back(data.index_inv[std::make_pair(i, j)]);
			while (id_layer > 0) {
				for (int k = 0;k < temp_edges.size();k++) {
					for (int m = 0;m < id_current_layer.size();m++) {
						if (temp_edges[k].second == id_current_layer[m]) {
							dependency_relationship[temp_edges[k].first][data.index_inv[std::make_pair(i, j)]] = true;
							temp_id_current_layer.push_back(temp_edges[k].first);
							break;
						}
					}
				}
				id_layer--;
				id_current_layer = temp_id_current_layer;
				temp_id_current_layer.clear();
			}
		}
	}
	cout << "the number of dependency edges: " << temp_edges.size() << endl;
	for (int i = 0;i < temp_collision_edges.size();i++) {
		if (dependency_relationship[temp_collision_edges[i].first][temp_collision_edges[i].second] == false) {
			dependency_relationship[temp_collision_edges[i].first][temp_collision_edges[i].second] = true;
			for (int j = 0;j < 10000;j++) {   //update
				if (dependency_relationship[temp_collision_edges[i].second][j] == true) {
					dependency_relationship[temp_collision_edges[i].first][j] = true;
				}
			}
			this->AddEdge_2(temp_collision_edges[i].first, temp_collision_edges[i].second);
			temp_edges.push_back(make_pair(temp_collision_edges[i].first, temp_collision_edges[i].second));
		}
	}
	cout << "the number of edges: " << temp_edges.size()<<endl;
	end_time = clock();
	cout << "&&&&&&& time of establish dependency graph (contain collision detection): " << double(end_time - start_time) / CLOCKS_PER_SEC << "s &&&&&&&" << endl;
}





void Layer_Graph::BuildDependencyGraph(std::vector<Eigen::Vector3d>& all_points)
{
	
	bool** Dependen_edges;
	Dependen_edges = new bool*[3000];
	for (int i = 0;i < 3000;i++)
		Dependen_edges[i] = new bool[3000];
	for (int i = 0;i < 3000;i++)
		for (int j = 0;j < 3000;j++)
			Dependen_edges[i][j] = false;
	std::vector<Eigen::Vector2d> all_2d_points;
	std::vector<Polygon> All_polygons;
 for (int i = 0;i < all_points.size();i++) {
		Eigen::Vector2d temp_point(all_points[i].x(), all_points[i].y());
		all_2d_points.push_back(temp_point);

		std::vector<Eigen::Vector2d> temp_2d_point;
		Eigen::Vector2d temp_2d_one_point = all_2d_points[i];
		double r = 0.2;
		//temp_2d_point.push_back(temp_2d_one_point);
		//temp_2d_one_point.x -= r;
		temp_2d_one_point.y() -= 2 * r;
		temp_2d_point.push_back(temp_2d_one_point);
		temp_2d_one_point.y() += 2 * r;
		temp_2d_point.push_back(temp_2d_one_point);
		temp_2d_one_point.y() += 2 * r;
		temp_2d_point.push_back(temp_2d_one_point);
		//temp_2d_one_point.x += r;
		//temp_2d_one_point.x += r;
		//temp_2d_one_point.y += r;
		
		/*for (double sita = 0;sita < 2 * PI;sita += 1) {
			temp_2d_one_point.x += r * sin(sita);
			temp_2d_one_point.y += r * cos(sita);
			temp_2d_point.push_back(temp_2d_one_point);
			temp_2d_one_point.x -= r * sin(sita);
			temp_2d_one_point.y -= r * cos(sita);
		}*/
		All_polygons.push_back(Polygon(ConstructPolygonPoints_2(temp_2d_point, dependence_offset/4*0.48)));
	}
	
	int cont_edges = 0;
 for (int i = 0;i < all_points.size();i++) {
		for (int j = 0;j < all_points.size();j++) {
			if (i != j && all_points[j].z() < all_points[i].z())
				if (All_polygons[j].JudgePointInside(all_2d_points[i])) {
					std::pair<int, int> Pair(i, j);
					Dependen_edges[i][j] = true;
					cont_edges++;
				}
		}
	}
	for (int i = 0;i < 3000;i++) {
		for (int j = 0;j < 3000;j++) {
			if (Dependen_edges[i][j] == true) {
				for (int t = 0;t < 3000;t++) {
					if (t != j && Dependen_edges[i][t] == true) {
						if (Dependen_edges[t][j] == true) {
							Dependen_edges[i][j] = false;
							break;
						}
					}
				}
			}
		}
	}

	ofstream all_edges("D:\\360MoveData\\Users\\zhong\\Desktop\\Others\\Open_Surface_Ceramics_Printing-master\\base_line\\all_edges.txt");
	all_edges << cont_edges << endl;;
	for (int i = 0;i < 3000;i++)
		for (int j = 0;j < 3000;j++)
			if (Dependen_edges[i][j] == true)
				all_edges << i << " " << j << endl;
	//Visual Vis;
	//Vis.generateModelForRendering_6(all_points, Dependen_edges);
}


void Layer_Graph::GetInitialOPP()
{
	clock_t start_time, end_time;
	start_time = clock();
	//std::cout << this->total_node_num << std::endl;
	int d_num = 0;
	for (int u = 0; u < this->total_node_num; u++) {
		for (int i = 0; i < G[u].size(); i++) {
			int v = this->edges[G[u][i]].GetTo(); d_num++;
			if (in_degree[v] >= 2 || out_degree[u] >= 2) {
				edge_deleted[G[u][i]] = true;
				d_num++;
			}
		}
	}
	std::cout << "d_num: " << d_num << std::endl;
	for (int i = 0; i < this->edges.size(); i++) {
		if (edge_deleted[i]) {
			int u = this->edges[i].GetFrom();
			int v = this->edges[i].GetTo();
			//std::cout << u  << " " << v << std::endl;
			in_degree[v]--;
			out_degree[u]--;
		}
	}
	for (int i = 0; i < this->total_node_num; i++) {
		if ((this->node_visited[i] == false) && (this->in_degree[i] == 0)) {
			//std::cout << i << std::endl;
			this->UpdateDegree(i, -1);
			std::vector<int> initial_opp;
			initial_opp.push_back(i);
			this->node_visited[i] = true;
			DFS(i, initial_opp);
			this->initial_opp_info.push_back(initial_opp);
		}
	}
	end_time = clock();
	std::cout << "***initial opp num***: " << this->initial_opp_info.size() << std::endl << std::endl;
	std::cout << "&&&&&&& time of establish initial opp graph: " << double(end_time - start_time) / CLOCKS_PER_SEC << "s &&&&&&&" << std::endl;
	//std::cout<<"&&&&&&& time of establish initial opp graph(vertex): " << double(end_time - start_time) / CLOCKS_PER_SEC << "s &&&&&&&" << endl;
#if DEFAULT_PRINTING_DIRECTION
	std::cout << "initial opp num: " << this->initial_opp_info.size() << std::endl;
#else

#endif
	this->OutputInitialOpp(file_name + "_initial_opp.txt");
}

bool Layer_Graph::IsDepend_collision(int i, int j)
{
	/////////*int id_top_layer_from = initial_opp_info[i][initial_opp_info[i].size() - 1];
	////////int id_bottom_layer_to = initial_opp_info[j][0];
	////////for (int k = cont_normal_dependency_edges;k < temp_edges.size();k++) {
	////////	if (temp_edges[k].first == id_top_layer_from && temp_edges[k].second == id_bottom_layer_to) {
	////////		return true;
	////////	}
	////////}
	////////return false;*/


	/*int layer_1, layer_2;
	for (int k = 0; k < data.slice_points[i].size(); k++) {
		layer_1 = data.index_inv[std::make_pair(i, k)];
		for (int l = 0; l < data.slice_points[j].size(); l++) {
			layer_2 = data.index_inv[std::make_pair(j, l)];
			for (int m = cont_normal_dependency_edges; m < temp_edges.size(); m++) {
				if (temp_edges[m].first == layer_1 && temp_edges[m].second == layer_2) {
					return true;
				}
			}
		}

	}*/

	
	for (int m = cont_normal_dependency_edges; m < temp_edges.size(); m++) {
		if (temp_edges[m].first == i && temp_edges[m].second == j) {
			return true;
		}
	}
	return false;
}

void Layer_Graph::DFS(int u, std::vector<int>& initial_opp)
{
	for (int i = 0; i < this->G[u].size(); i++) {
		int v = this->edges[G[u][i]].GetTo();
		if (edge_deleted[G[u][i]]) continue;
		if (this->node_visited[v]) continue;
		if (this->in_degree[v] != 0) continue;

		bool f1 = data.is_contour[data.index[u].first][data.index[u].second];
		bool f2 = data.is_contour[data.index[v].first][data.index[v].second];
		if (f1 != f2) {
			/*std::cout << data.index[u].first << " " << data.index[u].second << std::endl;
			std::cout << data.index[v].first << " " << data.index[v].second << std::endl;
			std::cout << *(data.slice_points[data.index[u].first][data.index[u].second].begin()) << std::endl;
			std::cout << *(data.slice_points[data.index[u].first][data.index[u].second].end() - 1) << std::endl;
			std::cout << *(data.slice_points[data.index[v].first][data.index[v].second].begin()) << std::endl;
			std::cout << *(data.slice_points[data.index[v].first][data.index[v].second].end() - 1) << std::endl;
			std::cout << "contour and segment " <<  f1 << " " << f2 << std::endl;*/
			continue;
		}
		this->UpdateDegree(v, -1);
		this->node_visited[v] = true;
		initial_opp.push_back(v);
		DFS(v, initial_opp);
		break;
	}
}

bool Layer_Graph::compare_two_node(std::vector<int>a, std::vector<int>b)
{
	bool jud_different = true;
	for (int k = 0; k < a.size();) {
		bool jud_different_2 = false;
		for (int l = 0; l < b.size();) {
			if (a[k] == b[l]) {
				k++;
				jud_different_2 = true;
				break;
			}
			else {
				l++;
			}
		}
		if (jud_different_2 == false) {
			jud_different = false;
			break;
		}
	}
	return jud_different;
}

void Layer_Graph::DFS_One(int u, std::vector<int>& medium_path, int num_blocks)
{
	for (int i = 0; i < this->total_node_num; i++) {
		int v = i;
		if (u == v) continue;
		if (this->node_visited[v]) continue;
		if (this->in_degree[v] != 0) continue;

		//bool jud_continue = false;
		//int sum_layers = 0;
		//for (int p = num_patches-1; p >= 0; p--) {    //�ⲿ����û��Ǹĳ�ÿ��patches֮�����������Ȼ�����д�
		//	sum_layers += cont_nodes_of_patches[p];
		//	if (i == sum_layers && medium_path.size() != 0) {
		//		jud_continue = true;
		//		for (int q = 0; q < medium_path.size(); q++) {
		//			if (medium_path[q] == sum_layers - 1)
		//				jud_continue = false;
		//		}
		//	}
		//}
		//if (jud_continue == true)
		//	continue;

		//if (i >= 166 && num_blocks <= 6)
			//continue;

		//////////collision dependency edges///////////////
		if (IsDepend_collision(u, v))
			continue;
		///////////////////////////////////////////////////

		bool jud_merge_layer = true;
		for (int k = 0; k < current_area_S.size(); k++) {
			if (has_subtractive_collision_dependency[current_area_S[k].id_layer][i] == true) {
				//update subtractive map
				Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
				for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
					if (temp_area_s_2.id_dep_layers[l] == i) {
						int index_id_dep_layers = l;
						for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
							int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
							int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
							for (int n = 0; n < current_point_map[k].size(); n++) {
								if (current_point_map[k][n].id_k == current_K) {
									current_point_map[k][n].flag_ori_unaccessible[current_ori] = true;
								}
							}
						}
					}
				}

				//judge whether can merge the layer
				bool jud_merge_layer_2 = true;
				for (int l = 0; l < current_point_map[k].size(); l++) {
					bool jud_merge_layer_3 = false;
					for (int m = 0; m < num_ori_sample; m++) {
						if (current_point_map[k][l].flag_ori_unaccessible[m] == false) {
							jud_merge_layer_3 = true;
						}
					}
					if (jud_merge_layer_3 == false)
						jud_merge_layer_2 = false;
				}
				if (jud_merge_layer_2 == false)
					jud_merge_layer = false;
			}
		}
		if (jud_merge_layer == true) {
			if (flag_layer_is_accessible[i] == false) {   //add area S to current block
				current_area_S.push_back(all_the_area_S[map_index_layers[i]]);

				////initial subtractive map
				vector<point_subtractive_map> temp_vec_point_subtractive_map;
				current_point_map.push_back(temp_vec_point_subtractive_map);
				for (int k = 0; k < all_the_area_S[map_index_layers[i]].all_k.size(); k++) {
					current_point_map[current_point_map.size() - 1].push_back(point_subtractive_map(all_the_area_S[map_index_layers[i]].all_k[k]));
				}

				for (int k = 0; k < all_the_area_S[map_index_layers[i]].id_dep_layers.size(); k++) {
					int id_layer_2 = all_the_area_S[map_index_layers[i]].id_dep_layers[k];
					has_subtractive_collision_dependency[i][id_layer_2] = true;
				}

				////add privious layers to current_point_map
				for (int k = 0; k < medium_path.size(); k++) {
					if (has_subtractive_collision_dependency[current_area_S[current_area_S.size() - 1].id_layer][medium_path[k]] == true) {
						Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[current_area_S.size() - 1].id_layer]];
						for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
							if (temp_area_s_2.id_dep_layers[l] == medium_path[k]) {
								int index_id_dep_layers = l;
								for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
									int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
									int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
									for (int n = 0; n < current_point_map[current_area_S.size() - 1].size(); n++) {
										if (current_point_map[current_area_S.size() - 1][n].id_k == current_K) {
											current_point_map[current_area_S.size() - 1][n].flag_ori_unaccessible[current_ori] = true;
										}
									}
								}
							}
						}
					}
				}
					
			}
			//all_blocks[num_blocks].push_back(i);
		}
		else {  //restore
			for (int k = 0; k < current_area_S.size(); k++) {
				if (has_subtractive_collision_dependency[current_area_S[k].id_layer][i] == true) {
					Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
					for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
						if (temp_area_s_2.id_dep_layers[l] == i) {
							int index_id_dep_layers = l;
							for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
								int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
								int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
								for (int n = 0; n < current_point_map[k].size(); n++) {
									if (current_point_map[k][n].id_k == current_K) {
										current_point_map[k][n].flag_ori_unaccessible[current_ori] = false;
									}
								}
							}
						}
					}

				}
			}

			continue;
		}

		this->node_visited[v] = true;
		this->UpdateDegree(v, -1);
		medium_path.push_back(v);

			
		DFS_One(v, medium_path,num_blocks);
		break;
	}
}


void Layer_Graph::DFS_ALL(int fa, int u, std::vector<int>& medium_path, int cont_num)
{
	//if (tree_node.size() >= 1000* cont_num)    //Ŀǰ�ܶ������
	//	return;


	bool go_on = false;
	for (int i = 0; i < this->total_node_num; i++) {
		int v = i;
		if (u == v) continue;
		if (this->node_visited[v]) continue;
		if (this->in_degree[v] != 0) continue;

		//////////collision dependency edges///////////////
		if (IsDepend_collision(u, v))
			continue;
		///////////////////////////////////////////////////

		bool jud_merge_layer = true;
		for (int k = 0; k < current_area_S.size(); k++) {
			if (has_subtractive_collision_dependency[current_area_S[k].id_layer][i] == true) {
				//update subtractive map
				Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
				for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
					if (temp_area_s_2.id_dep_layers[l] == i) {
						int index_id_dep_layers = l;
						for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
							int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
							int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
							for (int n = 0; n < current_point_map[k].size(); n++) {
								if (current_point_map[k][n].id_k == current_K) {
									current_point_map[k][n].flag_ori_unaccessible[current_ori] = true;
								}
							}
						}
					}
				}

				//judge whether can merge the layer
				bool jud_merge_layer_2 = true;
				for (int l = 0; l < current_point_map[k].size(); l++) {
					bool jud_merge_layer_3 = false;
					for (int m = 0; m < num_ori_sample; m++) {
						if (current_point_map[k][l].flag_ori_unaccessible[m] == false) {
							jud_merge_layer_3 = true;
						}
					}
					if (jud_merge_layer_3 == false)
						jud_merge_layer_2 = false;
				}
				if (jud_merge_layer_2 == false)
					jud_merge_layer = false;
			}
		}
		if (jud_merge_layer == true) {
			if (flag_layer_is_accessible[i] == false) {   //add area S to current block
				current_area_S.push_back(all_the_area_S[map_index_layers[i]]);

				////initial subtractive map
				vector<point_subtractive_map> temp_vec_point_subtractive_map;
				current_point_map.push_back(temp_vec_point_subtractive_map);
				for (int k = 0; k < all_the_area_S[map_index_layers[i]].all_k.size(); k++) {
					current_point_map[current_point_map.size() - 1].push_back(point_subtractive_map(all_the_area_S[map_index_layers[i]].all_k[k]));
				}

				for (int k = 0; k < all_the_area_S[map_index_layers[i]].id_dep_layers.size(); k++) {
					int id_layer_2 = all_the_area_S[map_index_layers[i]].id_dep_layers[k];
					has_subtractive_collision_dependency[i][id_layer_2] = true;
				}

				////add privious layers to current_point_map
				for (int k = 0; k < medium_path.size(); k++) {
					if (has_subtractive_collision_dependency[current_area_S[current_area_S.size() - 1].id_layer][medium_path[k]] == true) {
						Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[current_area_S.size() - 1].id_layer]];
						for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
							if (temp_area_s_2.id_dep_layers[l] == medium_path[k]) {
								int index_id_dep_layers = l;
								for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
									int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
									int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
									for (int n = 0; n < current_point_map[current_area_S.size() - 1].size(); n++) {
										if (current_point_map[current_area_S.size() - 1][n].id_k == current_K) {
											current_point_map[current_area_S.size() - 1][n].flag_ori_unaccessible[current_ori] = true;
										}
									}
								}
							}
						}
					}
				}

			}
			//all_blocks[num_blocks].push_back(i);
		}
		else {  //restore
			for (int k = 0; k < current_area_S.size(); k++) {
				if (has_subtractive_collision_dependency[current_area_S[k].id_layer][i] == true) {
					Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
					for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
						if (temp_area_s_2.id_dep_layers[l] == i) {
							int index_id_dep_layers = l;
							for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
								int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
								int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
								for (int n = 0; n < current_point_map[k].size(); n++) {
									if (current_point_map[k][n].id_k == current_K) {
										current_point_map[k][n].flag_ori_unaccessible[current_ori] = false;
									}
								}
							}
						}
					}

				}
			}
			continue;
		}

		this->node_visited[v] = true;
		this->UpdateDegree(v, -1);
		go_on = true;
		medium_path.push_back(v);
		DFS_ALL(fa, v, medium_path, cont_num);

		//restore
		medium_path.pop_back();
		for (int k = 0; k < current_area_S.size(); k++) {
			if (has_subtractive_collision_dependency[current_area_S[k].id_layer][i] == true) {
				Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
				for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
					if (temp_area_s_2.id_dep_layers[l] == i) {
						int index_id_dep_layers = l;
						for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
							int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
							int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
							for (int n = 0; n < current_point_map[k].size(); n++) {
								if (current_point_map[k][n].id_k == current_K) {
									current_point_map[k][n].flag_ori_unaccessible[current_ori] = false;
								}
							}
						}
					}
				}

			}
		}


		this->node_visited[v] = false;    
		this->UpdateDegree(v, 1);
	}
	if (go_on == false) {
		//ɾ������
		vector<int> temp_medium_path = medium_path;
		sort(temp_medium_path.begin(), temp_medium_path.end());
		bool jud_redundancy = false;
		for (int i = 0; i < temp_tree_node.size(); i++)
			if (temp_tree_node[i] == temp_medium_path)
				jud_redundancy = true;
		if (jud_redundancy == true)
			return;


		node++;
		tree_node.push_back(medium_path);
		temp_tree_node.push_back(temp_medium_path);
		pre_tree_index[node] = fa;
		Q_2.push(node);
	}
}

void Layer_Graph::BFS(std::vector<int>& medium_path, int num_blocks)
{
	int w = 10;
	int height = 0;
	std::vector<int> new_tree_node;
	int* new_pre_tree_index;
	new_pre_tree_index = new int[10000000];
	int cont_node = 0;
	vector<int> candidate_nodes;
	vector<int> terminate_nodes;
	memset(new_pre_tree_index, -1, sizeof(new_pre_tree_index));
	int root = -1;
	new_tree_node.push_back(root);
	node = 0;
	Q.push(node);
	while (Q.size() != 0 || candidate_nodes.size() != 0) {
		int now_node = Q.front();
		medium_path.clear();

		///////////////update/////////////
		while (new_pre_tree_index[now_node] != -1) {
			this->UpdateDegree(new_tree_node[now_node], -1);
			this->node_visited[new_tree_node[now_node]] = true;
			medium_path.push_back(new_tree_node[now_node]);
		
			if (new_tree_node.size() != 0) {
				for (int k = 0; k < current_area_S.size(); k++) {
					if (has_subtractive_collision_dependency[current_area_S[k].id_layer][new_tree_node[now_node]] == true) {
						//update subtractive map
						Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
						for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
							if (temp_area_s_2.id_dep_layers[l] == new_tree_node[now_node]) {
								int index_id_dep_layers = l;
								for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
									int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
									int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
									for (int n = 0; n < current_point_map[k].size(); n++) {
										if (current_point_map[k][n].id_k == current_K) {
											current_point_map[k][n].flag_ori_unaccessible[current_ori] = true;
										}
									}
								}
							}
						}
					}
				}
				if (flag_layer_is_accessible[new_tree_node[now_node]] == false) {   //add area S to current block
					current_area_S.push_back(all_the_area_S[map_index_layers[new_tree_node[now_node]]]);
					////initial subtractive map
					vector<point_subtractive_map> temp_vec_point_subtractive_map;
					current_point_map.push_back(temp_vec_point_subtractive_map);
					for (int k = 0; k < all_the_area_S[map_index_layers[new_tree_node[now_node]]].all_k.size(); k++) {
						current_point_map[current_point_map.size() - 1].push_back(point_subtractive_map(all_the_area_S[map_index_layers[new_tree_node[now_node]]].all_k[k]));
					}
					for (int k = 0; k < all_the_area_S[map_index_layers[new_tree_node[now_node]]].id_dep_layers.size(); k++) {
						int id_layer_2 = all_the_area_S[map_index_layers[new_tree_node[now_node]]].id_dep_layers[k];
						has_subtractive_collision_dependency[new_tree_node[now_node]][id_layer_2] = true;
					}
					////add privious layers to current_point_map
					for (int k = 0; k < medium_path.size(); k++) {
						if (has_subtractive_collision_dependency[current_area_S[current_area_S.size() - 1].id_layer][medium_path[k]] == true) {
							Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[current_area_S.size() - 1].id_layer]];
							for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
								if (temp_area_s_2.id_dep_layers[l] == medium_path[k]) {
									int index_id_dep_layers = l;
									for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
										int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
										int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
										for (int n = 0; n < current_point_map[current_area_S.size() - 1].size(); n++) {
											if (current_point_map[current_area_S.size() - 1][n].id_k == current_K) {
												current_point_map[current_area_S.size() - 1][n].flag_ori_unaccessible[current_ori] = true;
											}
										}
									}
								}
							}
						}
					}

				}
			}
			now_node = new_pre_tree_index[now_node];
		}
		//////////////////////////////////
		now_node = Q.front();
		bool jud_terminate = true;

		for (int i = 0; i < this->total_node_num; i++) {
			int v = i;
			if (new_tree_node.size()!=0 && new_tree_node[now_node] == v) continue;
			bool jud_continue_2 = false;
			for (int j = 0; j < candidate_nodes.size(); j++)
				if (new_tree_node[candidate_nodes[j]] == v) {
					jud_continue_2 = true;
					break;
				}
			if (jud_continue_2 == true) continue;
			if (this->node_visited[v]) continue;
			if (this->in_degree[v] != 0) continue;

			//////////collision dependency edges///////////////
			if (new_tree_node.size() != 0 && IsDepend_collision(new_tree_node[now_node], v))
				continue;
			///////////////////////////////////////////////////

			bool jud_merge_layer = true;
			for (int k = 0; k < current_area_S.size(); k++) {
				if (has_subtractive_collision_dependency[current_area_S[k].id_layer][i] == true) {
					//update subtractive map
					Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
					for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
						if (temp_area_s_2.id_dep_layers[l] == i) {
							int index_id_dep_layers = l;
							for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
								int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
								int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
								for (int n = 0; n < current_point_map[k].size(); n++) {
									if (current_point_map[k][n].id_k == current_K) {
										current_point_map[k][n].flag_ori_unaccessible[current_ori] = true;
									}
								}
							}
						}
					}

					//judge whether can merge the layer
					bool jud_merge_layer_2 = true;
					for (int l = 0; l < current_point_map[k].size(); l++) {
						bool jud_merge_layer_3 = false;
						for (int m = 0; m < num_ori_sample; m++) {
							if (current_point_map[k][l].flag_ori_unaccessible[m] == false) {
								jud_merge_layer_3 = true;
							}
						}
						if (jud_merge_layer_3 == false)
							jud_merge_layer_2 = false;
					}
					if (jud_merge_layer_2 == false)
						jud_merge_layer = false;
				}
			}
			if (jud_merge_layer == true) {
				if (flag_layer_is_accessible[i] == false) {   //add area S to current block
					current_area_S.push_back(all_the_area_S[map_index_layers[i]]);

					////initial subtractive map
					vector<point_subtractive_map> temp_vec_point_subtractive_map;
					current_point_map.push_back(temp_vec_point_subtractive_map);
					for (int k = 0; k < all_the_area_S[map_index_layers[i]].all_k.size(); k++) {
						current_point_map[current_point_map.size() - 1].push_back(point_subtractive_map(all_the_area_S[map_index_layers[i]].all_k[k]));
					}

					for (int k = 0; k < all_the_area_S[map_index_layers[i]].id_dep_layers.size(); k++) {
						int id_layer_2 = all_the_area_S[map_index_layers[i]].id_dep_layers[k];
						has_subtractive_collision_dependency[i][id_layer_2] = true;
					}

					////add privious layers to current_point_map
					for (int k = 0; k < medium_path.size(); k++) {
						if (has_subtractive_collision_dependency[current_area_S[current_area_S.size() - 1].id_layer][medium_path[k]] == true) {
							Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[current_area_S.size() - 1].id_layer]];
							for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
								if (temp_area_s_2.id_dep_layers[l] == medium_path[k]) {
									int index_id_dep_layers = l;
									for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
										int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
										int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
										for (int n = 0; n < current_point_map[current_area_S.size() - 1].size(); n++) {
											if (current_point_map[current_area_S.size() - 1][n].id_k == current_K) {
												current_point_map[current_area_S.size() - 1][n].flag_ori_unaccessible[current_ori] = true;
											}
										}
									}
								}
							}
						}
					}

				}
				//all_blocks[num_blocks].push_back(i);
			}
			else {  //restore
				for (int k = 0; k < current_area_S.size(); k++) {
					if (has_subtractive_collision_dependency[current_area_S[k].id_layer][i] == true) {
						Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
						for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
							if (temp_area_s_2.id_dep_layers[l] == i) {
								int index_id_dep_layers = l;
								for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
									int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
									int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
									for (int n = 0; n < current_point_map[k].size(); n++) {
										if (current_point_map[k][n].id_k == current_K) {
											current_point_map[k][n].flag_ori_unaccessible[current_ori] = false;
										}
									}
								}
							}
						}

					}
				}

				continue;
			}
			//this->node_visited[v] = false;  
			//this->UpdateDegree(v, 1);

			jud_terminate = false;
			node++;
			new_tree_node.push_back(v);
			new_pre_tree_index[node] = now_node;
			candidate_nodes.push_back(node);

			//restore
			for (int k = 0; k < current_area_S.size(); k++) {
				if (has_subtractive_collision_dependency[current_area_S[k].id_layer][i] == true) {
					Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
					for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
						if (temp_area_s_2.id_dep_layers[l] == i) {
							int index_id_dep_layers = l;
							for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
								int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
								int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
								for (int n = 0; n < current_point_map[k].size(); n++) {
									if (current_point_map[k][n].id_k == current_K) {
										current_point_map[k][n].flag_ori_unaccessible[current_ori] = false;
									}
								}
							}
						}
					}

				}
			}
		}

		if (jud_terminate == true)
			terminate_nodes.push_back(now_node);

		///////////////update/////////////
		while (new_pre_tree_index[now_node] != -1) {
			this->UpdateDegree(new_tree_node[now_node], 1);
			this->node_visited[new_tree_node[now_node]] = false;
			now_node = new_pre_tree_index[now_node];
		}

		//restore
		current_point_map.clear();
		current_area_S.clear();
		for (int j = 0; j < sum_layers; j++)
			for (int k = 0; k < sum_layers; k++)
				has_subtractive_collision_dependency[j][k] = false;
		//////////////////////////////////

		Q.pop();
		if (Q.size() == 0) {
			height++;
			int cont_w = 0;

			//sort    
			for (int i = 0; i < candidate_nodes.size(); i++)
				for (int j = i + 1; j < candidate_nodes.size(); j++) {
					if (flag_layer_is_accessible[new_tree_node[candidate_nodes[j]]] == false) {
						swap(candidate_nodes[i], candidate_nodes[j]);
						//swap(new_tree_node[i], new_tree_node[j]);
					}
						
				}

			while (candidate_nodes.size() != 0 && cont_w < w && cont_w < candidate_nodes.size()) {
				Q.push(candidate_nodes[cont_w]);
				cont_w++;
			}
			//if (cont_w >= w)
				//std::cout << "have reached W" << std::endl;
			//cout << "The number of candidates nodes:" << candidate_nodes.size() << endl;
			candidate_nodes.clear();
		}
	}


	//int now_node = node;    //��ѡһ��·��
	//std::vector<int> temp_1;
	//while (new_pre_tree_index[now_node] != -1) {
	//	temp_1.push_back(new_tree_node[now_node]);
	//	now_node = new_pre_tree_index[now_node];
	//}
	//std::reverse(temp_1.begin(), temp_1.end());
	//medium_path = temp_1;

	//for (int i = 0; i < temp_1.size(); i++) {
	//	this->node_visited[temp_1[i]] = true;
	//	this->UpdateDegree(temp_1[i], -1);
	//}


	vector<vector<int>> final_pathes;   //����·������ѡһ��
	vector<int> cont_area_S_for_pathes;
	for (int i = 0; i < terminate_nodes.size(); i++) {
		vector<int> current_path;
		current_path.clear();
		int cont_area_S_for_path = 0;
		int current_node = terminate_nodes[i];
		current_path.push_back(terminate_nodes[i]);
		if (flag_layer_is_accessible[new_tree_node[current_node]] == false)
			cont_area_S_for_path++;
		while (new_pre_tree_index[current_node] != 0) {
			current_node = new_pre_tree_index[current_node];
			current_path.push_back(current_node);
			if (flag_layer_is_accessible[new_tree_node[current_node]] == false)
				cont_area_S_for_path++;
		}
		std::reverse(current_path.begin(), current_path.end());
		final_pathes.push_back(current_path);
		cont_area_S_for_pathes.push_back(cont_area_S_for_path);
	}
	for (int i = 0; i < final_pathes.size(); i++) {
		for (int j = i + 1; j < final_pathes.size(); j++) {
			if (cont_area_S_for_pathes[j] >= cont_area_S_for_pathes[i]) {
				swap(cont_area_S_for_pathes[i], cont_area_S_for_pathes[j]);
				swap(final_pathes[i], final_pathes[j]);
			}
				
		}
	}
	vector<int> best_path;
	for (int i = 0; i < final_pathes[0].size(); i++) {
		best_path.push_back(new_tree_node[final_pathes[0][i]]);
	}
	medium_path = best_path;
	for (int i = 0; i < medium_path.size(); i++) {
		this->node_visited[medium_path[i]] = true;
		this->UpdateDegree(medium_path[i], -1);
	}
}


void Layer_Graph::OutputInitialOpp(const std::string& file_name)
{
	std::ofstream dstream(file_name.c_str());
	if (!dstream.is_open()) {
		std::cout << "can not open " << file_name << std::endl;
		return;
	}
	dstream << initial_opp_info.size() << std::endl;
	for (int i = 0; i < initial_opp_info.size(); i++) {
		dstream << initial_opp_info[i].size() << std::endl;
		int m, n;
		for (int j = 0; j < initial_opp_info[i].size(); j++) {
			m = data.index[initial_opp_info[i][j]].first;
			n = data.index[initial_opp_info[i][j]].second;
			dstream << data.slice_points[m][n].size() << std::endl;
			for (int k = 0; k < data.slice_points[m][n].size(); k++) {
                dstream << this->data.slice_points[m][n][k].x() << " " << this->data.slice_points[m][n][k].y() << " " << data.z_value[m][n][k] << std::endl;
			}
		}
	}
	dstream.close();
}

//void Layer_Graph::MappingBackLayers(vector<Eigen::Matrix3d> all_rotMatrix)
//{
//	int current_patch = num_patches-1;
//	int cont_layers = 0;
//	current_layers.clear();
//	current_layers.resize(data.slice_points.size());
//	
//	for (int i = 0; i < data.slice_points.size(); i++) {
//		if (cont_layers == cont_layers_of_patches[current_patch]) {
//			cont_layers = 0;
//			current_patch--;
//		}
//		current_layers[i].resize(data.slice_points[i].size());
//		for (int j = 0; j < data.slice_points[i].size(); j++) {
//			current_layers[i][j].resize(data.slice_points[i][j].size());
//			for (int k = 0; k < data.slice_points[i][j].size(); k++) {
//				current_layers[i][j][k].resize(3,1);
//				current_layers[i][j][k](0, 0) = data.slice_points[i][j][k].x;
//				current_layers[i][j][k](1, 0) = data.slice_points[i][j][k].y;
//				current_layers[i][j][k](2, 0) = data.z_value[i][j][k];
//				current_layers[i][j][k] = all_rotMatrix[num_patches - current_patch -1].inverse() * current_layers[i][j][k];
//			}
//		}	
//		cont_layers++;
//	}
//	//Visual Vis;
//	//Vis.generateModelForRendering(current_layers, file_name);
//}


void Layer_Graph::CollisionDetectionForAdditiveManufacturing(nozzle the_nozzle)
{
	clock_t start_time_9, end_time_9;
	
	vector<pair<int, int>> temp_collision_edges;
	temp_collision_edges.clear();
	
	//collision detect, add collision dependency edge
	for (int i = 0; i < data.slice_points.size(); i++) {
		for (int j = 0; j < data.slice_points[i].size(); j++) {
			for (int ii = i + 1; ii < data.slice_points.size(); ii++) {
				double circle_r;
				if ((ii - i) * dh < the_nozzle.nozzle_H_half)
					circle_r = the_nozzle.lowwer_surface_r + (ii - i) * dh * (the_nozzle.upper_surface_r - the_nozzle.lowwer_surface_r) / the_nozzle.nozzle_H_half;
				else
					circle_r = the_nozzle.upper_surface_r;

				for (int jj = 0; jj < data.slice_points[ii].size(); jj++) {
					bool jud_collision = false;
					if ((ii - i) * dh > the_nozzle.nozzle__H_total) //exceed nozzle_H
					{
						jud_collision = true;
						temp_collision_edges.push_back(make_pair(data.index_inv[std::make_pair(i, j)], data.index_inv[std::make_pair(ii, jj)]));
						continue;
					}
					jud_collision = ComponentBoundariesWithinDistance(
						data,
						i,
						j,
						ii,
						jj,
						circle_r * circle_r);
					if (jud_collision == true) {
						temp_collision_edges.push_back(make_pair(data.index_inv[std::make_pair(i, j)], data.index_inv[std::make_pair(ii, jj)]));
					}
				}
			}
		}
	}

	//start_time_9 = clock();
	//establish hase graph
	bool** dependency_relationship = new bool* [data.total_node_num];
	for (int i = 0; i < data.total_node_num; i++)
		dependency_relationship[i] = new bool[data.total_node_num];
	for (int i = 0; i < data.total_node_num; i++)
		for (int j = 0; j < data.total_node_num; j++)
			dependency_relationship[i][j] = false;
	//end_time_9 = clock();

	for (int i = 1; i < data.slice_points.size(); i++) {
		for (int j = 0; j < data.slice_points[i].size(); j++) {
			int id_layer = i;
			vector<int> id_current_layer;
			vector<int> temp_id_current_layer;
			id_current_layer.push_back(data.index_inv[std::make_pair(i, j)]);
			while (id_layer > 0) {
				for (int k = 0; k < temp_edges.size(); k++) {
					for (int m = 0; m < id_current_layer.size(); m++) {
						if (temp_edges[k].second == id_current_layer[m]) {
							dependency_relationship[temp_edges[k].first][data.index_inv[std::make_pair(i, j)]] = true;
							temp_id_current_layer.push_back(temp_edges[k].first);
							break;
						}
					}
				}
				id_layer--;
				id_current_layer = temp_id_current_layer;
				temp_id_current_layer.clear();
			}
		}
	}
	
	cont_normal_dependency_edges = temp_edges.size();

	
	//cout << "the number of dependency edges: " << temp_edges.size() << endl;
	for (int i = 0; i < temp_collision_edges.size(); i++) {
		if (dependency_relationship[temp_collision_edges[i].first][temp_collision_edges[i].second] == false) {
			dependency_relationship[temp_collision_edges[i].first][temp_collision_edges[i].second] = true;
			for (int j = 0; j < data.total_node_num; j++) {   //update
				if (dependency_relationship[temp_collision_edges[i].second][j] == true) {
					dependency_relationship[temp_collision_edges[i].first][j] = true;
				}
			}
			this->AddEdge_2(temp_collision_edges[i].first, temp_collision_edges[i].second);
			temp_edges.push_back(make_pair(temp_collision_edges[i].first, temp_collision_edges[i].second));
		}
	}
	
	for (int i = 0; i < data.total_node_num; i++)
		delete[]dependency_relationship[i];
	delete[]dependency_relationship;

	
	//cout << endl << "time..." << double(end_time_9 - start_time_9) / CLOCKS_PER_SEC << endl;
}

void merge_sort(vector<pair<int, int>> a, int l, int r) 
{
	if (l >= r) return;
	int mid = (l + r) >> 1;  
	merge_sort(a, l, mid);  
	merge_sort(a, mid + 1, r);
	int k = 0, i = l, j = mid + 1;
	vector<pair<int, int>> tmp;
	tmp.resize(a.size());
	while (i <= mid && j <= r) {
		if (a[i] < a[j]) tmp[k++] = a[i++];
		else tmp[k++] = a[j++];
	}
	while (i <= mid) tmp[k++] = a[i++];
	while (j <= r) tmp[k++] = a[j++];


	for (i = l, j = 0; i <= r; i++, j++)	a[i] = tmp[j];
}



void Layer_Graph::CollisionDetectionForSubtractiveManufacturing(nozzle the_nozzle)
{
	cout << "Doing collision detection for subtractive manufacturing......" << endl;

	sum_layers = 0;
	int step = 10;   //10
	string file_name_2 = file_name;
	int cont_num = 0;
	SAMPLE_ON_BALL sampling;
	sampling.OrientationSamplePoints();
	Eigen::Matrix3d rotMatrix;
	std::vector<std::vector<std::vector<bool>>> flag_accessible_points;
	vector<pair<int, int>> temp_collision_edges;
	vector<int> id_orientation;
	vector<int> id_k;
	double nozzle_par = (the_nozzle.upper_surface_r - the_nozzle.lowwer_surface_r) / the_nozzle.nozzle_H_half;

	flag_accessible_points.resize(current_layers.size());
	for (int i = 0; i < current_layers.size(); i++) {
		flag_accessible_points[i].resize(current_layers[i].size());
		for (int j = 0; j < current_layers[i].size(); j++) {
			flag_accessible_points[i][j].resize(current_layers[i][j].size());
			sum_layers++;
			for (int k = 0; k < current_layers[i][j].size(); k++)
				flag_accessible_points[i][j][k] = false;
		}
	}

	for (int ori = 0; ori < sampling.sample_points.size(); ori++) {
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling.sample_points[ori]);
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		std::vector<std::vector<std::vector<Eigen::MatrixXd>>> temp_layers = current_layers;

		for (int i = 0; i < temp_layers.size(); i++)
			for (int j = 0; j < temp_layers[i].size(); j++)
				for (int k = 0; k < temp_layers[i][j].size(); k++)
					temp_layers[i][j][k] = rotMatrix.inverse() * temp_layers[i][j][k];
		
		//file_name_2 = file_name + to_string(cont_num);



		///////////////////collision detection////////////////////////
		int cont_accessible_points = 0, cont_unaccessible_points = 0;
		double circle_r;

		for (int i = 0; i < temp_layers.size(); i++) {     
			for (int j = 0; j < temp_layers[i].size(); j++) {
				for (int k = 0; k < temp_layers[i][j].size(); k += step) {   //�²���
					if (flag_accessible_points[i][j][k] == true) {
						cont_accessible_points++;
						continue;
					}
					cont_unaccessible_points++;
					bool jud_collision = false;
					for (int ii = 0; ii < temp_layers.size(); ii++) {
						for (int jj = 0; jj < temp_layers[ii].size(); jj++) {
							if (i != ii || j != jj) {
								bool jud_collision_2 = false;
								for (int kk = 0; kk < temp_layers[ii][jj].size(); kk += step) {
									//if (i != ii || j != jj || k != kk) {
									double diff = temp_layers[ii][jj][kk](2, 0) - temp_layers[i][j][k](2, 0);
									if (diff <= 0)
										continue;
									else if (diff < the_nozzle.nozzle_H_half)
										circle_r = the_nozzle.lowwer_surface_r + (temp_layers[ii][jj][kk](2, 0) - temp_layers[i][j][k](2, 0)) * nozzle_par;
									else
										circle_r = the_nozzle.upper_surface_r;

									if (diff > the_nozzle.nozzle__H_total) {  //exceed nozzle_H
										jud_collision_2 = true;
										break;
									}
									else if (pow(temp_layers[ii][jj][kk](0, 0) - temp_layers[i][j][k](0, 0), 2) + pow(temp_layers[ii][jj][kk](1, 0) - temp_layers[i][j][k](1, 0), 2) - pow(circle_r, 2) < 0) {
										jud_collision_2 = true;
										break;
									}
									//}
								}
								/*if (jud_collision == true)
									break;*/
								if (jud_collision_2 == true) {
									jud_collision = true;
									//temp_collision_edges.push_back(make_pair(data.index_inv[std::make_pair(i, j)], data.index_inv[std::make_pair(ii, jj)]));
									//id_orientation.push_back(ori);
									//id_k.push_back(k);
								}
							}
						}
					}
					if (jud_collision == false) {
						flag_accessible_points[i][j][k] = true;
					}
				}
			}
		}
		cout << "id of orientation:" << ori << endl;
		cout << "number of accessible points:" << cont_accessible_points << endl;
		cout << "number of unaccessible points:" << cont_unaccessible_points << endl <<endl;
						//break;
	}


	//*******************************************//
	/////////////////find area S//////////////////
	//*******************************************//
	for (int ori = 0; ori < sampling.sample_points.size(); ori++) {
		Eigen::Vector3d vectorBefore(0, 0, 1);
		Eigen::Vector3d vectorAfter(sampling.sample_points[ori]);
		rotMatrix = Eigen::Quaterniond::FromTwoVectors(vectorBefore, vectorAfter).toRotationMatrix();
		std::vector<std::vector<std::vector<Eigen::MatrixXd>>> temp_layers = current_layers;

		for (int i = 0; i < temp_layers.size(); i++)
			for (int j = 0; j < temp_layers[i].size(); j++)
				for (int k = 0; k < temp_layers[i][j].size(); k++)
					temp_layers[i][j][k] = rotMatrix.inverse() * temp_layers[i][j][k];

		ofstream vis_unaccessible_points(file_name+"_unaccessible_points.obj");
		file_name_2 = file_name + to_string(cont_num);

		std::vector<Eigen::MatrixXd> vis_points;
		double circle_r;
		for (int i = 0; i < temp_layers.size(); i++) {
			for (int j = 0; j < temp_layers[i].size(); j++) {
				for (int k = 0; k < temp_layers[i][j].size(); k += step) {   //�²���
					if (flag_accessible_points[i][j][k] == true) {
						continue;
					}
					vis_unaccessible_points << "v " << current_layers[i][j][k](0,0) << " " << current_layers[i][j][k](1, 0) << " " << current_layers[i][j][k](2, 0) << endl;
					vis_points.push_back(current_layers[i][j][k]);

					
					for (int ii = 0; ii < temp_layers.size(); ii++) {
						for (int jj = 0; jj < temp_layers[ii].size(); jj++) {
							if (i != ii || j != jj) {
								bool jud_collision_2 = false;   //���λ��Ӧ�ö���
								for (int kk = 0; kk < temp_layers[ii][jj].size(); kk += step) {
									//if (i != ii || j != jj || k != kk) {
									double diff = temp_layers[ii][jj][kk](2, 0) - temp_layers[i][j][k](2, 0);
									if (diff <= 0)
										continue;
									else if (diff < the_nozzle.nozzle_H_half)
										circle_r = the_nozzle.lowwer_surface_r + (temp_layers[ii][jj][kk](2, 0) - temp_layers[i][j][k](2, 0)) * nozzle_par;
									else
										circle_r = the_nozzle.upper_surface_r;

									if (diff > the_nozzle.nozzle__H_total) {  //exceed nozzle_H
										jud_collision_2 = true;
										break;
									}
									else if (pow(temp_layers[ii][jj][kk](0, 0) - temp_layers[i][j][k](0, 0), 2) + pow(temp_layers[ii][jj][kk](1, 0) - temp_layers[i][j][k](1, 0), 2) - pow(circle_r, 2) < 0) {
										jud_collision_2 = true;
										break;
									}
									//}
								}
								/*if (jud_collision == true)
									break;*/
								if (jud_collision_2 == true) {
									temp_collision_edges.push_back(make_pair(data.index_inv[std::make_pair(i, j)], data.index_inv[std::make_pair(ii, jj)]));
									id_orientation.push_back(ori);
									id_k.push_back(k);
								}
							}
						}
					}
				}
			}
		}

		Visual Vis;
		//Vis.generateModelForRendering(temp_layers, file_name_2);
		//Vis.generateModelForRendering_3(vectorAfter, file_name_2, vis_points);
		//break;
		if(ori == 0)
			creat_ball(file_name,vis_points);
		cont_num++;
	}
	//*******************************************//


	temp_area_S the_temp_area_S(temp_collision_edges.size());
	for (int i = 0; i < temp_collision_edges.size(); i++) {
		the_temp_area_S.id_layers[i] = temp_collision_edges[i].first;
		the_temp_area_S.id_points[i] = id_k[i];
		the_temp_area_S.id_ori[i] = id_orientation[i];
		the_temp_area_S.subtractive_collision_edges[i] = temp_collision_edges[i];
	}



	////////////************generate Subtractive collision dependency edges*********////////////////////////
	for (int i = 0; i < temp_collision_edges.size(); i++) {
		bool jud_break = false;
		for (int j = 0; j < all_the_area_S.size(); j++) {
			if (the_temp_area_S.id_layers[i] == all_the_area_S[j].id_layer) {
				jud_break = true;
				bool jud_break_2 = false;
				for (int k = 0; k < all_the_area_S[j].id_dep_layers.size(); k++) {
					if (the_temp_area_S.subtractive_collision_edges[i].second == all_the_area_S[j].id_dep_layers[k]) {
						jud_break_2 = true;
						all_the_area_S[j].K_and_ori[k].push_back(make_pair(the_temp_area_S.id_points[i], the_temp_area_S.id_ori[i]));
						break;
					}
				}
				if (jud_break_2 == false) {
					int temp_id_dep_layers = the_temp_area_S.subtractive_collision_edges[i].second;
					vector<pair<int, int>> temp_K_and_ori;
					all_the_area_S[j].id_dep_layers.push_back(temp_id_dep_layers);
					all_the_area_S[j].K_and_ori.push_back(temp_K_and_ori);
					j--;
				}
				else
					break;
			}
		}
		if (jud_break == false) {
			Area_S temp_area;
			temp_area.id_layer = the_temp_area_S.id_layers[i];
			all_the_area_S.push_back(temp_area);
			i--;
		}
	}

	for (int i = 0; i < all_the_area_S.size(); i++) {
		for(int j =0;j<all_the_area_S[i].K_and_ori.size();j++)
			for (int k = 0; k < all_the_area_S[i].K_and_ori[j].size(); k++) {
				bool jud_exist = false;
				for (int l = 0; l < all_the_area_S[i].all_k.size(); l++) {
					if (all_the_area_S[i].all_k[l] == all_the_area_S[i].K_and_ori[j][k].first) {
						jud_exist = true;
					}
				}
				if (jud_exist == false) {
					all_the_area_S[i].all_k.push_back(all_the_area_S[i].K_and_ori[j][k].first);
				}
			}
	}
	cout << "Collision Detection For Subtractive Manufacturing over" <<endl;
}

void Layer_Graph::MergeLayersWithDefaultOrder()
{ 
	flag_layer_is_accessible = new bool[sum_layers];
	for (int i = 0; i < sum_layers; i++)
		flag_layer_is_accessible[i] = true;
	for (int i = 0; i < all_the_area_S.size(); i++) {
		flag_layer_is_accessible[all_the_area_S[i].id_layer] = false;
		map_index_layers.insert({ all_the_area_S[i].id_layer, i });
	}
	

	///////**************default searching order**************/////////
	///��ʱ����������ײ������
	///�������ɺϲ�layerʱ���ֿ�
	
	has_subtractive_collision_dependency = new bool* [sum_layers];
	for (int i = 0; i < sum_layers; i++)
		has_subtractive_collision_dependency[i] = new bool[sum_layers];
	for (int i = 0; i < sum_layers; i++)
		for (int j = 0; j < sum_layers; j++)
			has_subtractive_collision_dependency[i][j] = false;

	vector<int> temp_vec_all_blocks;
	all_blocks.push_back(temp_vec_all_blocks);
	int num_blocks = 0;
	for(int i =0;i<current_layers.size();i++)
		for (int j = 0; j < current_layers[i].size(); j++) {    //Ĭ��˳�������Ƿ���֧��������ϵ
			int id_layer = data.index_inv[std::make_pair(i, j)];
			bool jud_merge_layer = true;
			for (int k = 0; k < current_area_S.size(); k++) {
				if (has_subtractive_collision_dependency[current_area_S[k].id_layer][id_layer] == true) {
					//update subtractive map
					Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[k].id_layer]];
					for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
						if (temp_area_s_2.id_dep_layers[l] == id_layer) {
							int index_id_dep_layers = l;
							for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
								int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
								int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
								for (int n = 0; n < current_point_map[k].size(); n++) {
									if (current_point_map[k][n].id_k == current_K) {
										current_point_map[k][n].flag_ori_unaccessible[current_ori] = true;
									}
								}
							}
						}
					}

					//judge whether can merge the layer
					bool jud_merge_layer_2 = true;
					for (int l = 0; l < current_point_map[k].size(); l++) {
						bool jud_merge_layer_3 = false;
						for (int m = 0; m < num_ori_sample; m++) {
							if (current_point_map[k][l].flag_ori_unaccessible[m] == false) {
								jud_merge_layer_3 = true;
							}
						}
						if (jud_merge_layer_3 == false)
							jud_merge_layer_2 = false;
					}
					if (jud_merge_layer_2 == false)
						jud_merge_layer = false;
				}
			}
			if (jud_merge_layer == true) {
				if (flag_layer_is_accessible[id_layer] == false) {   //add area S to current block
					current_area_S.push_back(all_the_area_S[map_index_layers[id_layer]]);	

					////initial subtractive map
					vector<point_subtractive_map> temp_vec_point_subtractive_map;
					current_point_map.push_back(temp_vec_point_subtractive_map);
					for (int k = 0; k < all_the_area_S[map_index_layers[id_layer]].all_k.size(); k++) {
						current_point_map[current_point_map.size() - 1].push_back(point_subtractive_map(all_the_area_S[map_index_layers[id_layer]].all_k[k]));
					}

					for (int k = 0; k < all_the_area_S[map_index_layers[id_layer]].id_dep_layers.size(); k++) {
						int id_layer_2 = all_the_area_S[map_index_layers[id_layer]].id_dep_layers[k];
						has_subtractive_collision_dependency[id_layer][id_layer_2] = true;
					}

					////add privious layers to current_point_map
					for (int k = 0; k < all_blocks[num_blocks].size(); k++) {
						if (has_subtractive_collision_dependency[current_area_S[current_area_S.size() - 1].id_layer][all_blocks[num_blocks][k]] == true) {
							Area_S temp_area_s_2 = all_the_area_S[map_index_layers[current_area_S[current_area_S.size() - 1].id_layer]];
							for (int l = 0; l < temp_area_s_2.id_dep_layers.size(); l++) {
								if (temp_area_s_2.id_dep_layers[l] == all_blocks[num_blocks][k]) {
									int index_id_dep_layers = l;
									for (int m = 0; m < temp_area_s_2.K_and_ori[index_id_dep_layers].size(); m++) {
										int current_K = temp_area_s_2.K_and_ori[index_id_dep_layers][m].first;
										int current_ori = temp_area_s_2.K_and_ori[index_id_dep_layers][m].second;
										for (int n = 0; n < current_point_map[current_area_S.size() - 1].size(); n++) {
											if (current_point_map[current_area_S.size() - 1][n].id_k == current_K) {
												current_point_map[current_area_S.size() - 1][n].flag_ori_unaccessible[current_ori] = true;
											}
										}
									}
								}
							}
						}
					}
				}
				all_blocks[num_blocks].push_back(id_layer);
			}

			else if (jud_merge_layer == false) {    //creat new block
				////clear subtractive map
				current_point_map.clear();
				current_area_S.clear();

				for (int i = 0; i < sum_layers; i++)
					for (int j = 0; j < sum_layers; j++)
						has_subtractive_collision_dependency[i][j] = false;

				num_blocks++;
				vector<int> temp_vec_all_blocks;
				all_blocks.push_back(temp_vec_all_blocks);
				j--;
			}
		}
	cout << "**************************************" << endl;
	cout << "Number of blocks:" << all_blocks.size() << endl;;
	cout << "**************************************" << endl;

	int cont_all_blocks = 0;
	for (int i = 0; i < all_blocks.size(); i++)
		for (int j = 0; j < all_blocks[i].size(); j++)
			cont_all_blocks++;

	cout << cont_all_blocks << " " << sum_layers << endl;

	std::vector<std::vector<Eigen::MatrixXd>> all_lines;
	for (int i = 0; i < current_layers.size(); i++) {
		for (int j = 0; j < current_layers[i].size(); j++) {
			std::vector<Eigen::MatrixXd> temp_lines;
			all_lines.push_back(temp_lines);
			for (int k = 0; k < current_layers[i][j].size(); k++) {
				all_lines[all_lines.size() - 1].push_back(current_layers[i][j][k]);
			}
		}
	}

	Visual Vis;
	Vis.generateModelForRendering_4(all_lines, file_name, all_blocks);
}

void Layer_Graph::MergeLayersWithGreedyAlgorithm(int* cont_layers_of_patches)
{
	/*int sum_layers_2 = 0;
	cont_nodes_of_patches = new int[num_patches];
	for (int i = num_patches - 1; i >= 0; i--)
		cont_nodes_of_patches[i] = 0;
	for (int i = num_patches - 1; i >= 0; i--) {
		for (int j = sum_layers_2; j < sum_layers_2+cont_layers_of_patches[i]; j++) {
			for (int k = 0; k < current_layers[j].size(); k++) {
				cont_nodes_of_patches[i]++;
			}
		}
		sum_layers_2 += cont_layers_of_patches[i];
	}*/


	flag_layer_is_accessible = new bool[sum_layers];
	for (int i = 0; i < sum_layers; i++)
		flag_layer_is_accessible[i] = true;
	for (int i = 0; i < all_the_area_S.size(); i++) {
		flag_layer_is_accessible[all_the_area_S[i].id_layer] = false;
		map_index_layers.insert({ all_the_area_S[i].id_layer, i });
	}

	has_subtractive_collision_dependency = new bool* [sum_layers];
	for (int i = 0; i < sum_layers; i++)
		has_subtractive_collision_dependency[i] = new bool[sum_layers];
	for (int i = 0; i < sum_layers; i++)
		for (int j = 0; j < sum_layers; j++)
			has_subtractive_collision_dependency[i][j] = false;

	//vector<int> temp_vec_all_blocks_2;
	//all_blocks.push_back(temp_vec_all_blocks_2);
	int num_blocks = 0;
	
////////////////////////////////////////////**************greedy search**************////////////////////////////////////////////////
	for (int i = 0; i < this->total_node_num; i++) {
		if (this->node_visited[i] == false && this->in_degree[i] == 0) {
			this->UpdateDegree(i, -1);
			vector<int> temp_vec_all_blocks;
			temp_vec_all_blocks.push_back(i);
			this->node_visited[i] = true;
			DFS_One(i, temp_vec_all_blocks,num_blocks);
			all_blocks.push_back(temp_vec_all_blocks);
			//all_blocks.push_back(temp_vec_all_blocks);
			cout << "num_blocks" << num_blocks << ",size:  " << temp_vec_all_blocks.size() << endl;
			num_blocks++;

			//restore
			current_point_map.clear();
			current_area_S.clear();
			for (int j = 0; j < sum_layers; j++)
				for (int k = 0; k < sum_layers; k++)
					has_subtractive_collision_dependency[j][k] = false;
		}
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "**************************************" << endl;
	cout << "Number of blocks:" << all_blocks.size() << endl;;
	cout << "**************************************" << endl;

	int cont_all_blocks = 0;
	for (int i = 0; i < all_blocks.size(); i++)
		for (int j = 0; j < all_blocks[i].size(); j++)
			cont_all_blocks++;

	cout << cont_all_blocks << " " << sum_layers << endl;

	std::vector<std::vector<Eigen::MatrixXd>> all_lines;
	for (int i = 0; i < current_layers.size(); i++) {
		for (int j = 0; j < current_layers[i].size(); j++) {
			std::vector<Eigen::MatrixXd> temp_lines;
			all_lines.push_back(temp_lines);
			for (int k = 0; k < current_layers[i][j].size(); k++) {
				all_lines[all_lines.size() - 1].push_back(current_layers[i][j][k]);
			}
		}
	}

	Visual Vis;
	Vis.generateModelForRendering_4(all_lines, file_name, all_blocks);
}


void Layer_Graph::MergeLayersWithBeamSearch(int* cont_layers_of_patches)
{
	flag_layer_is_accessible = new bool[sum_layers];
	for (int i = 0; i < sum_layers; i++)
		flag_layer_is_accessible[i] = true;
	for (int i = 0; i < all_the_area_S.size(); i++) {
		flag_layer_is_accessible[all_the_area_S[i].id_layer] = false;
		map_index_layers.insert({ all_the_area_S[i].id_layer, i });
	}

	has_subtractive_collision_dependency = new bool* [sum_layers];
	for (int i = 0; i < sum_layers; i++)
		has_subtractive_collision_dependency[i] = new bool[sum_layers];
	for (int i = 0; i < sum_layers; i++)
		for (int j = 0; j < sum_layers; j++)
			has_subtractive_collision_dependency[i][j] = false;

	int num_blocks = 0;
	///////////////////////////////////////////////beam search//////////////////////////////////////////////////
	clock_t start_time, end_time;
	start_time = clock();
	int w = 1;
	int height = 0;
	pre_tree_index = new int[10000000];
	tree_node.clear();
	temp_tree_node.clear();
	memset(pre_tree_index, -1, sizeof(pre_tree_index));
	std::vector<int> root;
	tree_node.push_back(root);
	temp_tree_node.push_back(root);
	node = 0;
	Q.push(node);
	int cont_num = 1;
	while (Q.size() != 0 || Q_2.size() != 0) {
		int now_node = Q.front();
		while (pre_tree_index[now_node] != -1) {
			for (int i = 0; i < tree_node[now_node].size(); i++) {
				this->UpdateDegree(tree_node[now_node][i], -1);
				this->node_visited[tree_node[now_node][i]] = true;
			}
			now_node = pre_tree_index[now_node];
		}
		now_node = Q.front();
		for (int i = 0; i < this->total_node_num; i++) {
			if (this->node_visited[i] == false && this->in_degree[i] == 0) {
				std::vector<int> temp_vec_all_blocks;
				temp_vec_all_blocks.push_back(i);
				this->UpdateDegree(i, -1);
				this->node_visited[i] = true;
				this->DFS_ALL(now_node, i, temp_vec_all_blocks, cont_num);
				cont_num++;

						//for (int j = 1; j < temp_vec_all_blocks.size(); j++) {   //����ɾ���������ӣ�
						//	this->node_visited[temp_vec_all_blocks[j]] = false;
						//	this->UpdateDegree(temp_vec_all_blocks[j], 1);
						//}

				//restore
				this->UpdateDegree(i, 1);
				this->node_visited[i] = false;
				current_point_map.clear();
				current_area_S.clear();
				for (int j = 0; j < sum_layers; j++)
					for (int k = 0; k < sum_layers; k++)
						has_subtractive_collision_dependency[j][k] = false;
			}
		}

		while (pre_tree_index[now_node] != -1) {
			for (int i = 0; i < tree_node[now_node].size(); i++) {
				this->UpdateDegree(tree_node[now_node][i], 1);
				this->node_visited[tree_node[now_node][i]] = false;
			}
			now_node = pre_tree_index[now_node];
		}
		Q.pop();
		if (Q.size() == 0) {
			height++;
			int cont_w = 0;
			while (Q_2.size() != 0 && cont_w < w) {
				Q.push(Q_2.front());
				Q_2.pop();
				cont_w++;
			}
			if (cont_w >= w)
				std::cout << "have reached W" << std::endl;
			Q_2 = std::queue<int>();
		}
	}

	int min_path = MAX_I;
	int cnt1, cnt2;
	int now_node;
	for (int i = node; i >= 0; i--) {
		now_node = i;
		cnt1 = cnt2 = 0;
		while (pre_tree_index[now_node] != -1) {
			cnt1++;
			for (int j = 0; j < tree_node[now_node].size(); j++) {
				cnt2++;
			}
			now_node = pre_tree_index[now_node];
		}
		if (cnt2 == this->total_node_num) min_path = std::min(cnt1, min_path);
	}
	std::cout << "first path cover medium opp number_search_tree: " << min_path << std::endl;
	for (int i = node; i >= 0; i--) {
		now_node = i;
		cnt2 = 0;
		std::vector<std::vector<int>> temp_1;
		while (pre_tree_index[now_node] != -1) {
			std::vector<int> temp_2;
			for (int j = 0; j < tree_node[now_node].size(); j++) {
				temp_2.push_back(tree_node[now_node][j]);
				cnt2++;
			}
			temp_1.push_back(temp_2);
			now_node = pre_tree_index[now_node];
		}
		std::reverse(temp_1.begin(), temp_1.end());
		if (cnt2 == this->total_node_num && temp_1.size() == min_path) {
			this->all_solutions.push_back(temp_1);
		}
	}

	for (int i = 0; i < all_solutions.size(); i++)
		for (int j = i + 1; j < all_solutions.size(); j++) {
			bool jud_different = true;
			for (int k = 0; k < all_solutions[i].size();) {
				bool jud_different_2 = false;
				for (int l = 0; l < all_solutions[j].size();) {
					if (compare_two_node(all_solutions[i][k], all_solutions[j][l]) == true) {
						k++;
						jud_different_2 = true;
						break;
					}
					else {
						l++;
					}
				}
				if (jud_different_2 == false) {
					jud_different = false;
					break;
				}

			}
			if (jud_different == true) {
				std::vector<std::vector<std::vector<int>>>::iterator itor;
				int cont_itor = 0;
				for (itor = all_solutions.begin(); itor != all_solutions.end(); itor++) {
					if (cont_itor == i) {
						all_solutions.erase(itor);
						break;
					}
					cont_itor++;
				}
				i--;
				//cont_needless_node++;
				break;
			}
		}

	end_time = clock();
	std::cout << "&&&&&&& time of PC-MPC calculation: " << double(end_time - start_time) / CLOCKS_PER_SEC << "s &&&&&&&" << std::endl;

	delete[] pre_tree_index;
	all_blocks = all_solutions[0];
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}


void Layer_Graph::MergeLayersWithHeuristicSearch(int* cont_layers_of_patches)
{

	flag_layer_is_accessible = new bool[sum_layers];
	for (int i = 0; i < sum_layers; i++)
		flag_layer_is_accessible[i] = true;
	for (int i = 0; i < all_the_area_S.size(); i++) {
		flag_layer_is_accessible[all_the_area_S[i].id_layer] = false;
		map_index_layers.insert({ all_the_area_S[i].id_layer, i });
	}

	has_subtractive_collision_dependency = new bool* [sum_layers];
	for (int i = 0; i < sum_layers; i++)
		has_subtractive_collision_dependency[i] = new bool[sum_layers];
	for (int i = 0; i < sum_layers; i++)
		for (int j = 0; j < sum_layers; j++)
			has_subtractive_collision_dependency[i][j] = false;

	//vector<int> temp_vec_all_blocks_2;
	//all_blocks.push_back(temp_vec_all_blocks_2);
	int num_blocks = 0;

	////////////////////////////////////////////**************Heuristic Search**************////////////////////////////////////////////////
	for (int i = 0; i < this->total_node_num; i++) {
		if (this->node_visited[i] == false && this->in_degree[i] == 0) {
			//this->UpdateDegree(i, -1);
			vector<int> temp_vec_all_blocks;
			//temp_vec_all_blocks.push_back(i);
			//this->node_visited[i] = true;
			BFS(temp_vec_all_blocks, num_blocks);
			all_blocks.push_back(temp_vec_all_blocks);
			//all_blocks.push_back(temp_vec_all_blocks);
			cout << "num_blocks" << num_blocks << ",size:  " << temp_vec_all_blocks.size() << endl;
			num_blocks++;

			//restore
			current_point_map.clear();
			current_area_S.clear();
			for (int j = 0; j < sum_layers; j++)
				for (int k = 0; k < sum_layers; k++)
					has_subtractive_collision_dependency[j][k] = false;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "**************************************" << endl;
	cout << "Number of blocks:" << all_blocks.size() << endl;;
	cout << "**************************************" << endl;

	int cont_all_blocks = 0;

	for (int i = 0; i < all_blocks.size(); i++)
		for (int j = 0; j < all_blocks[i].size(); j++)
			cont_all_blocks++;

	cout << cont_all_blocks << " " << sum_layers << endl;

	std::vector<std::vector<Eigen::MatrixXd>> all_lines;
	for (int i = 0; i < current_layers.size(); i++) {
		for (int j = 0; j < current_layers[i].size(); j++) {
			std::vector<Eigen::MatrixXd> temp_lines;
			all_lines.push_back(temp_lines);
			for (int k = 0; k < current_layers[i][j].size(); k++) {
				all_lines[all_lines.size() - 1].push_back(current_layers[i][j][k]);
			}
		}
	}

	Visual Vis;
	Vis.generateModelForRendering_4(all_lines, file_name, all_blocks);
}
