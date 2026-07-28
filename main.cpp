#include <filesystem>
//#include <igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include<igl/readSTL.h>
#include<igl/writeOBJ.h>
#include <cstdio>
#include <cstdlib>
#include "stl2obj/vectornd.h"
#include "stl2obj/geometry.h"
#include "stl2obj/importstl.h"
#include "stl2obj/exportobj.h"
#include <Eigen/Dense>
#include <tuple>
#include<iostream>
#include<igl/writeSTL.h>
#include "helpers.h"
#include "data.h"
#include "sample_on_ball.h"
#include "HybridManufacturing.h"
#include <ranges>
#include <vector>
#include <io.h>
#include<direct.h> 

#include <CGAL/Circular_kernel_intersections.h>
#include <CGAL/Exact_circular_kernel_2.h>

namespace fs = std::filesystem;

enum class SubtractiveDecompositionMode
{
	BeamSearch = 0,
	Global = 1,
	Local = 2,
	LocalInitializedGlobal = 3,
	RandomMultiStartGlobal = 4,
	RandomMultiStartLocal = 5
};

#ifndef VASCO_ROOT_DIR
#define VASCO_ROOT_DIR "."
#endif


// ���ߺ�������
void get_need_file(string path, vector<string>& file, vector<string>& file_name, string ext);
void change_name(string my_file_path, string suff);
vector<fs::path> collect_files(string dir);
void remove_files(vector<fs::path> files, string dir);
nozzle create_nozzle();
cutter create_cutter(float tolerance);
bool process_model(
	const std::string& my_file_path,
	const std::string& suff,
	const nozzle& the_nozzle,
	const cutter& cutting_tool,
	SubtractiveDecompositionMode subtractive_decomposition_mode,
	int subtractive_decomposition_height,
	int subtractive_random_restart_count,
	unsigned int subtractive_random_seed,
	int subtractive_random_thread_count,
	double subtractive_block_boundary_overlap_cost_per_unit);


//clock_t start_time, end_time;

int main()
{
	std::cout << "[Debug] Current Working Directory: "
		<< std::filesystem::current_path() << std::endl;

	const std::string iniPath = "input_files.ini";
	const IniData ini = LoadIni(iniPath);
	const std::string my_file_path = GetIniString(ini, "Input", "BaseDir", "data\\transformed_mesh_10K");
	const std::vector<std::string> my_file_name = GetIniStringList(ini, "Input", "FileNames");
	int subtractive_decomposition_mode_value =
		GetIniInt(ini, "Run", "SubtractiveDecompositionMode", 0);
	if (subtractive_decomposition_mode_value < static_cast<int>(SubtractiveDecompositionMode::BeamSearch)
		|| subtractive_decomposition_mode_value > static_cast<int>(SubtractiveDecompositionMode::RandomMultiStartLocal)) {
		std::cerr << "[main] Invalid SubtractiveDecompositionMode="
			<< subtractive_decomposition_mode_value
			<< "; expected 0 (beam search), 1 (global), 2 (local),"
			<< " 3 (local-initialized global), 4 (random-multistart global),"
			<< " or 5 (random-multistart local). Use 0."
			<< std::endl;
		subtractive_decomposition_mode_value =
			static_cast<int>(SubtractiveDecompositionMode::BeamSearch);
	}
	const auto subtractive_decomposition_mode =
		static_cast<SubtractiveDecompositionMode>(subtractive_decomposition_mode_value);
	const int subtractive_decomposition_height =
		GetIniInt(ini, "Run", "SubtractiveDecompositionHeight", 3);
	const int subtractive_random_restart_count =
		std::max(1, GetIniInt(ini, "Run", "SubtractiveRandomRestartCount", 10));
	const int subtractive_random_seed_value =
		std::max(0, GetIniInt(ini, "Run", "SubtractiveRandomSeed", 12345));
	const unsigned int subtractive_random_seed =
		static_cast<unsigned int>(subtractive_random_seed_value);
	const int subtractive_random_thread_count =
		std::max(0, GetIniInt(ini, "Run", "SubtractiveRandomThreadCount", 0));
	const double subtractive_block_boundary_overlap_cost_per_unit =
		std::max(
			0.0,
			GetIniDouble(
				ini,
				"Run",
				"SubtractiveBlockBoundaryOverlapCostPerUnit",
				10.0));

	if (my_file_name.empty()) {
		std::cerr << "[main] No input file names found in " << iniPath << std::endl;
		return 1;
	}

	const float tolerance = 0.05f;
	const nozzle the_nozzle = create_nozzle();
	const cutter cutting_tool = create_cutter(tolerance);

	for (const auto& suff : my_file_name) {
		process_model(
			my_file_path,
			suff,
			the_nozzle,
			cutting_tool,
			subtractive_decomposition_mode,
			subtractive_decomposition_height,
			subtractive_random_restart_count,
			subtractive_random_seed,
			subtractive_random_thread_count,
			subtractive_block_boundary_overlap_cost_per_unit);
	}

	return 0;
}

// �����ļ�����ָ���ļ���׺���ļ�
void get_need_file(string path, vector<string>& file, vector<string>& file_name, string ext)
{

	fs::path dir(path);
	std::error_code ec;

	// Ŀ¼У��
	if (!fs::exists(dir, ec) || !fs::is_directory(dir, ec)) {
		std::cerr << "[get_need_file] Directory not found or not a directory: " << path
			<< " (" << ec.message() << ")" << std::endl;
		return;
	}

	fs::directory_iterator it(dir, ec);
	if (ec) {
		std::cerr << "[get_need_file] Iterate failed: " << path
			<< " (" << ec.message() << ")" << std::endl;
		return;
	}

	std::vector<std::pair<std::string, std::string>> tmp; // {fullpath, filename}
	for (const auto& entry : it) {
		if (!entry.is_regular_file(ec)) continue;
		const auto& p = entry.path();

		std::string e = p.extension().string();

		if (e == ext) {
			tmp.emplace_back(p.string(), p.filename().string());
		}
	}

	file.clear();
	file_name.clear();

	file.reserve(tmp.size());
	file_name.reserve(tmp.size());
	for (auto& pr : tmp) {
		file.emplace_back(std::move(pr.first));
		file_name.emplace_back(std::move(pr.second));
	}
}
// ��0_0
void change_name(string my_file_path, string suff)
{
	vector<string> my_file;
	vector<string> my_file_name;

	get_need_file(my_file_path, my_file, my_file_name, suff);
	for (int i = 0; i < my_file.size(); i++)
	{
		string pre_file_name = my_file_name[i].substr(0, my_file_name[i].find("."));


		string newname = my_file_path + '\\' + pre_file_name + "-0_0" + suff;
		// int rename(char *oldname, char *newname);
		rename(my_file[i].c_str(), newname.c_str());

	}
}
// �ҵ�Ŀ���ļ����е������ļ�
vector<fs::path> collect_files(string dir)
{
	auto files = vector<fs::path>();
	auto iterator = fs::directory_iterator(dir);

	for (auto p : iterator)
	{
		auto file = p.path();
		if (file.extension() == ".obj" || file.extension() == ".stl" || file.extension() == ".gcode" || file.extension() == ".txt")
			files.push_back(p.path());
	}

	return files;
}
// �ƶ��ļ�
void remove_files(vector<fs::path> files, string dir)
{
	fs::path new_dir = fs::path(dir);
	if (!fs::exists(dir))
	{
		fs::create_directory(new_dir);
	}

	for (auto p : files)
	{
		fs::path new_path = new_dir / p.filename();
		fs::rename(p, new_path);
	}
}




nozzle create_nozzle()
{
	nozzle the_nozzle;
	the_nozzle.upper_surface_r = 4.5;
	the_nozzle.lowwer_surface_r = 1.0;
	the_nozzle.nozzle__H_total = 15;
	the_nozzle.nozzle_H_half = 5;
	return the_nozzle;
}

cutter create_cutter(float tolerance)
{
	cutter cutting_tool;
	cutting_tool.cylinder_r = 1.5;
	cutting_tool.cylinder_height = 25;
	cutting_tool.ball_r = 1.5;
	cutting_tool.carriage_r = 20;
	cutting_tool.carriage_height = 30;

	cutting_tool.cylinder_r *= (1 + tolerance);
	cutting_tool.carriage_r *= (1 + tolerance);
	cutting_tool.cylinder_height *= (1 - tolerance);
	cutting_tool.carriage_height *= (1 - tolerance);
	cutting_tool.ball_r *= (1 + tolerance);
	return cutting_tool;
}

bool process_model(
	const std::string& my_file_path,
	const std::string& suff,
	const nozzle& the_nozzle,
	const cutter& cutting_tool,
	SubtractiveDecompositionMode subtractive_decomposition_mode,
	int subtractive_decomposition_height,
	int subtractive_random_restart_count,
	unsigned int subtractive_random_seed,
	int subtractive_random_thread_count,
	double subtractive_block_boundary_overlap_cost_per_unit)
{
	std::string file_name = my_file_path + "\\" + suff + "\\" + suff;
	std::string path_obj = file_name + ".obj";
	std::cout << path_obj << std::endl;

	ifstream file_normal(file_name + ".txt");
	int num_patches = 0;
	file_normal >> num_patches;
	vector<Eigen::Matrix3d> all_rotMatrix;
	(void)num_patches;
	(void)all_rotMatrix;

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd N;

	igl::readOBJ(path_obj, V, F);
	vasco::core::getMeshNormal(V, F, N);

	Katana::Instance().stl.saveStlFromObj(file_name + "-0_0" + ".stl", V, F);
	igl::writeOBJ(file_name + "-0_0" + ".obj", V, F); //�¼ӵ���䣬�����Ժ�beamsearchȱ��obj������
	HybridManufacturing hybrid_manufacturing(my_file_path, suff, V, F, N);
	hybrid_manufacturing.graph_cut_block_boundary_overlap_cost_per_unit =
		subtractive_block_boundary_overlap_cost_per_unit;
	hybrid_manufacturing.open_vis_voronoi = 0;
	hybrid_manufacturing.open_vis_red_points = 1;
	hybrid_manufacturing.open_vis_green_points = 0;
	hybrid_manufacturing.open_vis_stair_case = 1;
	hybrid_manufacturing.open_vis_additive_accessibility_debug = 0;
	hybrid_manufacturing.open_vis_additive_self_support_debug = 1;
	hybrid_manufacturing.open_change_orientation = 0;

	if (subtractive_decomposition_mode == SubtractiveDecompositionMode::Global) {
		cutter cutting_tool_sub = cutting_tool;
		std::cout << "[main] SubtractiveDecompositionMode=1 (global), height="
			<< subtractive_decomposition_height << std::endl;
		hybrid_manufacturing.subtractive_accessibility_decomposition_global(
			subtractive_decomposition_height,
			cutting_tool_sub);
		return true;
	}
	if (subtractive_decomposition_mode == SubtractiveDecompositionMode::Local) {
		cutter cutting_tool_sub = cutting_tool;
		std::cout << "[main] SubtractiveDecompositionMode=2 (local), height="
			<< subtractive_decomposition_height << std::endl;
		const int local_label_count =
			hybrid_manufacturing.subtractive_accessibility_decomposition_local(
				subtractive_decomposition_height,
				cutting_tool_sub);
		std::cout << "[main] Local subtractive decomposition label count="
			<< local_label_count << std::endl;
		return local_label_count != 999;
	}
	if (subtractive_decomposition_mode == SubtractiveDecompositionMode::LocalInitializedGlobal) {
		cutter cutting_tool_sub = cutting_tool;
		std::cout << "[main] SubtractiveDecompositionMode=3 (local-initialized global), height="
			<< subtractive_decomposition_height << std::endl;
		hybrid_manufacturing.subtractive_accessibility_decomposition_local_initialized_global(
			subtractive_decomposition_height,
			cutting_tool_sub);
		return true;
	}
	if (subtractive_decomposition_mode == SubtractiveDecompositionMode::RandomMultiStartGlobal) {
		cutter cutting_tool_sub = cutting_tool;
		std::cout << "[main] SubtractiveDecompositionMode=4 (random-multistart global), height="
			<< subtractive_decomposition_height
			<< ", restarts=" << subtractive_random_restart_count
			<< ", seed=" << subtractive_random_seed
			<< ", threads=" << subtractive_random_thread_count
			<< std::endl;
		hybrid_manufacturing.subtractive_accessibility_decomposition_random_multistart_global(
			subtractive_decomposition_height,
			cutting_tool_sub,
			subtractive_random_restart_count,
			subtractive_random_seed,
			subtractive_random_thread_count);
		return true;
	}
	if (subtractive_decomposition_mode == SubtractiveDecompositionMode::RandomMultiStartLocal) {
		cutter cutting_tool_sub = cutting_tool;
		std::cout << "[main] SubtractiveDecompositionMode=5 (random-multistart local), height="
			<< subtractive_decomposition_height
			<< ", restarts=" << subtractive_random_restart_count
			<< ", seed=" << subtractive_random_seed
			<< ", threads=" << subtractive_random_thread_count
			<< std::endl;
		const int local_label_count =
			hybrid_manufacturing.subtractive_accessibility_decomposition_random_multistart_local(
				subtractive_decomposition_height,
				cutting_tool_sub,
				subtractive_random_restart_count,
				subtractive_random_seed,
				subtractive_random_thread_count);
		std::cout << "[main] Random-multistart local subtractive decomposition label count="
			<< local_label_count << std::endl;
		return local_label_count != 999;
	}

	hybrid_manufacturing.InitializeVoronoi();
	int flag = hybrid_manufacturing.CollisionDetectionForSubtractiveManufacturing(cutting_tool);

	//std::string src_file = my_file_path + "\\" + suff;
	//std::string target_file1 = "data\\done\\" + suff;
	//fs::create_directory(target_file1);
	//if (flag == -1) {
	//	auto files = collect_files(src_file);
	//	remove_files(files, target_file1);
	//	fs::remove(src_file);
	//	return false;
	//}

	hybrid_manufacturing.outer_beam_search(the_nozzle, cutting_tool);
	return true;
}
