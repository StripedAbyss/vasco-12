for (int block_index = 1; block_index <= height_of_beam_search; ++block_index) {


	Slicer_2 slicer_load;
	slicer_load.load((file_name + "-" + to_string(block_index - 1) + "_0" + ".obj").c_str());
	const Slicer_2& slicer = slicer_load;	//加个const看看有没有涉及更改
	if (block_index == 1)
		slicer_load_current_patch.load(".\\vis\\block_patch-" + to_string(block_index) + "_" + ".obj");
	else
		slicer_load_current_patch = slicer_load_next_patch;
	
	if (block_index + 1 <= height_of_beam_search)
		slicer_load_next_patch.load(".\\vis\\block_patch-" + to_string(block_index + 1) + "_" + ".obj");
	else {
		slicer_load_next_patch.positions.clear();
		slicer_load_next_patch.triangles.clear();
	}

	//// Register the mesh with Polyscope
	//if (slicer_load_current_patch.triangles.size() > 0)
	//	polyscope::registerSurfaceMesh("current patch", slicer_load_current_patch.positions, slicer_load_current_patch.triangles);

	//// Register the mesh with Polyscope
	//if (slicer_load_next_patch.triangles.size() > 0)
	//	polyscope::registerSurfaceMesh("next patch", slicer_load_next_patch.positions, slicer_load_next_patch.triangles);

	// Show the GUI
	//polyscope::show();

	if (slicer_load_current_patch.triangles.size() == 0)
		continue;
	/////////////////Get sample points//////////////////
	vector<vasco::core::Vec3> sample_points_current(slicer_load_current_patch.triangles.size());
	for (int i = 0; i < slicer_load_current_patch.triangles.size(); i++) {
		auto i_1 = slicer_load_current_patch.triangles[i][0];
		auto i_2 = slicer_load_current_patch.triangles[i][1];
		auto i_3 = slicer_load_current_patch.triangles[i][2];

		auto V_1 = slicer_load_current_patch.positions[i_1];
		auto V_2 = slicer_load_current_patch.positions[i_2];
		auto V_3 = slicer_load_current_patch.positions[i_3];

		double a = distance3d(V_1, V_2);
		double b = distance3d(V_1, V_3);
		double c = distance3d(V_2, V_3);

		vasco::core::Vec3 V_incentre;

		if (a + b + c == 0) {
			V_incentre = V_1;
		}
		else
			for (int j = 0; j < 3; ++j) {
				V_incentre[j] = (a * V_1[j] + b * V_2[j] + c * V_3[j]) / (a + b + c);
			}
		sample_points_current[i] = V_incentre;
	}

	vector<vasco::core::Vec3> sample_points_next(slicer_load_next_patch.triangles.size());
	for (int i = 0; i < slicer_load_next_patch.triangles.size(); i++) {
		auto i_1 = slicer_load_next_patch.triangles[i][0];
		auto i_2 = slicer_load_next_patch.triangles[i][1];
		auto i_3 = slicer_load_next_patch.triangles[i][2];

		auto V_1 = slicer_load_next_patch.positions[i_1];
		auto V_2 = slicer_load_next_patch.positions[i_2];
		auto V_3 = slicer_load_next_patch.positions[i_3];

		double a = distance3d(V_1, V_2);
		double b = distance3d(V_1, V_3);
		double c = distance3d(V_2, V_3);

		vasco::core::Vec3 V_incentre;

		if (a + b + c == 0) {
			V_incentre = V_1;
		}
		else
			for (int j = 0; j < 3; ++j) {
				V_incentre[j] = (a * V_1[j] + b * V_2[j] + c * V_3[j]) / (a + b + c);
			}
		sample_points_next[i] = V_incentre;
	}


	cout << "% " << slicer_load.triangles.size() << endl;
	cout << "% " << sample_points_current.size() << endl;
	cout << "% " << sample_points_next.size() << endl;

	vector<vector<int>> accessible_ori_current = getAccessOri(slicer, slicer_load_current_patch, sample_points_current, cutting_tool);
	vector<vector<int>> accessible_ori_next = getAccessOri(slicer, slicer_load_next_patch, sample_points_next, cutting_tool);

	vector<vector<Eigen::MatrixXd>> vis_points;
	vector<vector<vector<Eigen::Vector3d>>> vis_lines;

	//////////////////////////////////graph cut////////////////////////////////////////
	//不应该出现cont_ori == 0的情况
	int cont_revise = 0;
	for (int i = 0; i < accessible_ori_current.size(); i++) {
		int cont_ori = 0;
		for (int j = 0; j < accessible_ori_current[i].size(); j++) {
			if (accessible_ori_current[i][j] == 0)
				cont_ori++;
		}
		//暂时先强制任意方向可达
		if (cont_ori == 0) {
			cont_revise++;
			//cout << "*** " << need_detect_triangle[i][0] << endl;
			for (int j = 0; j < accessible_ori_current[i].size(); j++)
				accessible_ori_current[i][j] = 0;
		}
		/*if (accessible_ori_of_need_detect_V[i][123] != 0)
			cout << "no" << endl;*/
	}
	cout << cont_revise << endl;

	vector<int> accessible_index_in_next_patch;
	//把next patch里面不可达的剃掉，不参与graph cut
	cont_revise = 0;
	for (int i = 0; i < accessible_ori_next.size(); i++) {
		int cont_ori = 0;
		for (int j = 0; j < accessible_ori_next[i].size(); j++) {
			if (accessible_ori_next[i][j] == 0)
				cont_ori++;
		}
		if (cont_ori > 0) {
			accessible_index_in_next_patch.push_back(i);
		}
		else {
			++cont_revise;
		}
		/*if (accessible_ori_of_need_detect_V[i][123] != 0)
			cout << "no" << endl;*/
	}
	cout << cont_revise << endl;

	vector<vector<int>> accessible_ori_of_need_detect_V;
	for (const auto& it : accessible_ori_current) {
		accessible_ori_of_need_detect_V.push_back(it);
	}
	for (const auto& index : accessible_index_in_next_patch) {
		accessible_ori_of_need_detect_V.push_back(accessible_ori_next[index]);
	}

	vector<std::pair<int,int>> pixels_relations;
	vector<int> length_edges;
	for (int i = 0; i < slicer_load_current_patch.triangles.size(); i++)
		for (int j = 0; j < slicer_load_current_patch.triangles.size(); j++) {
			bool jud_adjacent = false;
			vector<int> temp_edge(2);
			for (int ii = 0; ii < 3; ii++) {
				for (int jj = 0; jj < 3; jj++) {
					if (slicer_load_current_patch.triangles[i][ii] == slicer_load_current_patch.triangles[j][jj]) {
						/*temp_edge[0] = i;
						temp_edge[1] = j;
						pixels_relations.push_back(temp_edge);
						break;*/
						if (slicer_load_current_patch.triangles[i][(ii + 1) % 3] == slicer_load_current_patch.triangles[j][(jj + 1) % 3] ||
							slicer_load_current_patch.triangles[i][(ii + 1) % 3] == slicer_load_current_patch.triangles[j][(jj + 2) % 3]) {
							pixels_relations.push_back({ i,j });
							int index1 = slicer_load_current_patch.triangles[i][ii];
							int index2 = slicer_load_current_patch.triangles[i][(ii + 1) % 3];
							float temp_length = distanceVec3(slicer_load_current_patch.positions[index1], slicer.positions[index2]) * 100;
							length_edges.push_back(temp_length);
							jud_adjacent = true;
							break;
						}
						else if (slicer_load_current_patch.triangles[i][(ii + 2) % 3] == slicer_load_current_patch.triangles[j][(jj + 1) % 3] ||
							slicer_load_current_patch.triangles[i][(ii + 2) % 3] == slicer_load_current_patch.triangles[j][(jj + 2) % 3]) {
							temp_edge[0] = i;
							temp_edge[1] = j;
							pixels_relations.push_back({ i,j });
							int index1 = slicer_load_current_patch.triangles[i][ii];
							int index2 = slicer_load_current_patch.triangles[i][(ii + 2) % 3];
							float temp_length = distanceVec3(slicer_load_current_patch.positions[index1], slicer.positions[index2]) * 100;
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
		}

	for (int i = 0; i < accessible_index_in_next_patch.size(); i++)
		for (int j = 0; j < accessible_index_in_next_patch.size(); j++) {
			vasco::core::Tri3 tri_i = slicer_load_next_patch.triangles[accessible_index_in_next_patch[i]];
			vasco::core::Tri3 tri_j = slicer_load_next_patch.triangles[accessible_index_in_next_patch[j]];
			bool jud_adjacent = false;
			vector<int> temp_edge(2);
			for (int ii = 0; ii < 3; ii++) {
				for (int jj = 0; jj < 3; jj++) {
					if (tri_i[ii] == tri_j[jj]) {
						/*temp_edge[0] = i;
						temp_edge[1] = j;
						pixels_relations.push_back(temp_edge);
						break;*/
						if (tri_i[(ii + 1) % 3] == tri_j[(jj + 1) % 3] ||
							tri_i[(ii + 1) % 3] == tri_j[(jj + 2) % 3]) {
							int vertex_i = accessible_index_in_next_patch[i] + slicer_load_current_patch.triangles.size();
							int vertex_j = accessible_index_in_next_patch[j] + slicer_load_current_patch.triangles.size();
							pixels_relations.push_back({ vertex_i , vertex_j });
							int index1 = tri_i[ii];
							int index2 = tri_i[(ii + 1) % 3];
							float temp_length = distanceVec3(slicer_load_next_patch.positions[index1], slicer_load_next_patch.positions[index2]) * 10;
							length_edges.push_back(temp_length);
							jud_adjacent = true;
							break;
						}
						else if (tri_i[(ii + 2) % 3] == tri_j[(jj + 1) % 3] ||
							tri_i[(ii + 2) % 3] == tri_j[(jj + 2) % 3]) {
							int vertex_i = accessible_index_in_next_patch[i] + slicer_load_current_patch.triangles.size();
							int vertex_j = accessible_index_in_next_patch[j] + slicer_load_current_patch.triangles.size();
							pixels_relations.push_back({ vertex_i , vertex_j });
							int index1 = tri_i[ii];
							int index2 = tri_i[(ii + 2) % 3];
							float temp_length = distanceVec3(slicer_load_next_patch.positions[index1], slicer_load_next_patch.positions[index2]) * 10;
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
		}


	int sum_of_2_size = accessible_ori_current.size() + accessible_index_in_next_patch.size();
	cout << endl << "Graph Cut For 2................" << block_index << endl;
	cout << "% " << sum_of_2_size << endl;
	cout << "% " << sampling_subtractive.sample_points.size() << endl;


	vector<int> result = GeneralGraph_DArraySArraySpatVarying2(sum_of_2_size, sampling_subtractive.sample_points.size(), accessible_ori_of_need_detect_V, pixels_relations, length_edges);

	vector<int> all_result_labels = result;
	std::sort(all_result_labels.begin(), all_result_labels.end());
	auto temp_it = std::unique(all_result_labels.begin(), all_result_labels.end());
	all_result_labels.resize(std::distance(all_result_labels.begin(), temp_it));


	vector<vector<vasco::core::Tri3>> vis_triangles;
	vector<vasco::core::Vec3> vis_positions;

	for (const auto& pos : slicer_load_current_patch.positions) {
		vis_positions.push_back({ pos[0], pos[1] + 140.0, pos[2] });
	}

	for (const auto& pos : slicer_load_next_patch.positions) {
		vis_positions.push_back({ pos[0], pos[1] + 140.0, pos[2] });
	}

	auto current_pos_cnt = slicer_load_current_patch.positions.size();

	int cont_patch = 0;
	vector<Eigen::Vector3d> points_in_cell;
	vector<cv::Point3d> normals;
	//cout << "normal:" << result[0]<<endl;
	ofstream ofile(".\\vis\\normal_of_pathch2es" + std::to_string(block_index) + ".txt");
	for (int label : all_result_labels) {

		vector<Eigen::MatrixXd> temp_vecc(0);
		vis_points.push_back(temp_vecc);
		vector<vector<Eigen::Vector3d>> temp_vec(0);
		vis_lines.push_back(temp_vec);
		vector<vasco::core::Tri3> temp_tri(0);
		vis_triangles.push_back(temp_tri);
		normals.push_back(sampling_subtractive.sample_points[label]);
		ofile << normals[normals.size() - 1].x << " " << normals[normals.size() - 1].y << " " << normals[normals.size() - 1].z << endl;

		Eigen::Vector3d first_point(0.0, 0.0, 0.0);
		for (int ijj = 0; ijj < result.size(); ++ijj) {
			if (result[ijj] != label) {
				continue;
			}

			int index = ijj;
			if (ijj < accessible_ori_current.size()) {

				if (accessible_ori_current[index][label] != 0) {
					std::cout << "Error: accessible_ori_current[index][label] != 0 at index " << index << " and label " << label << std::endl;
				}


				first_point = Eigen::Vector3d(slicer_load_current_patch.positions[slicer_load_current_patch.triangles[index][0]][0],
					slicer_load_current_patch.positions[slicer_load_current_patch.triangles[index][0]][1],
					slicer_load_current_patch.positions[slicer_load_current_patch.triangles[index][0]][2]);

				vector<Eigen::Vector3d> temp_vec;

				temp_vec.clear();
				//cout << index_V_need_detect[i] << endl;
				for (int j = 0; j < 3; j++) {
					Eigen::Vector3d temp_vec_2(slicer_load_current_patch.positions[slicer_load_current_patch.triangles[index][j]][0],
						slicer_load_current_patch.positions[slicer_load_current_patch.triangles[index][j]][1],
						slicer_load_current_patch.positions[slicer_load_current_patch.triangles[index][j]][2]);
					temp_vec.push_back(temp_vec_2);
				}
				vis_lines[cont_patch].push_back(temp_vec);
				vis_triangles.rbegin()->push_back(slicer_load_current_patch.triangles[index]);

			}
			else {
				index -= int(accessible_ori_current.size());
				index = accessible_index_in_next_patch[index];

				if (index < 0 || index >= slicer_load_next_patch.triangles.size()) {
					std::cout << "Index out of bounds: " << index << " for slicer_load_next_patch.triangles of size " << slicer_load_next_patch.triangles.size() << std::endl;
				}

				if (accessible_ori_next[index][label] != 0) {
					std::cout << "Error: accessible_ori_next[index][label] != 0 at index " << index << " and label " << label << std::endl;
				}


				double p1 = slicer_load_next_patch.positions[slicer_load_next_patch.triangles[index][0]][0];
				double p2 = slicer_load_next_patch.positions[slicer_load_next_patch.triangles[index][0]][1];
				double p3 = slicer_load_next_patch.positions[slicer_load_next_patch.triangles[index][0]][2];
				first_point = Eigen::Vector3d(p1, p2, p3);

				vector<Eigen::Vector3d> temp_vec;

				temp_vec.clear();
				//cout << index_V_need_detect[i] << endl;
				for (int j = 0; j < 3; j++) {
					Eigen::Vector3d temp_vec_2(slicer_load_next_patch.positions[slicer_load_next_patch.triangles[index][j]][0],
						slicer_load_next_patch.positions[slicer_load_next_patch.triangles[index][j]][1],
						slicer_load_next_patch.positions[slicer_load_next_patch.triangles[index][j]][2]);
					temp_vec.push_back(temp_vec_2);
				}
				vis_lines[cont_patch].push_back(temp_vec);
				auto tri_index = slicer_load_next_patch.triangles[index];
				tri_index[0] += current_pos_cnt;
				tri_index[1] += current_pos_cnt;
				tri_index[2] += current_pos_cnt;

				vis_triangles.rbegin()->push_back(tri_index);
			}



		}
		points_in_cell.push_back(first_point);
		cont_patch++;
	}

	if (accessible_index_in_next_patch.size() != 0)
	{
		size_t next_face_cnt = slicer_load_next_patch.triangles.size();
		vector<vasco::core::Tri3> new_triangles;
		auto iter = accessible_index_in_next_patch.begin();
		for (size_t i = 0; i < next_face_cnt; ++i) {
			if (iter != accessible_index_in_next_patch.end() && *iter != i) {
				new_triangles.push_back(slicer_load_next_patch.triangles[i]);
			}
			if (iter != accessible_index_in_next_patch.end() && *iter == i) {
				++iter;
			}
		}
		slicer_load_next_patch.triangles = new_triangles;
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
		std::ofstream dstream(".\\vis\\2patchtwo-" + to_string(block_index) + "_" + to_string(t) + ".stl");
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
	ofile.close();


	for (int t = 0; t < vis_triangles.size(); t++) {
		std::string mesh_name = "2subtractive_2patch_" + to_string(block_index) + "_" + to_string(t);
		polyscope::registerSurfaceMesh(mesh_name, vis_positions, vis_triangles[t]);
	}
	polyscope::show();


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