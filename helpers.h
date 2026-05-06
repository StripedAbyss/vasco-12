#pragma once
#ifndef HELPERS_H_
#define HELPERS_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <filesystem>

#include <Eigen/Dense>
using namespace std;


struct nozzle
{
	double upper_surface_r;
	double lowwer_surface_r;
	double nozzle__H_total;
	double nozzle_H_half;
};

struct cutter
{
	double cylinder_r;
	double cylinder_height;
	double ball_r;
	double carriage_r;
	double carriage_height;
	double cylinder_height_threshold;
	double carriage_check_radius_sq;
	double cylinder_check_radius_sq;
	double cylinder_r_sq;
	double carriage_r_sq;
	double total_height;
};

const int terminate_threshold_of_number_of_faces = 1000;  //1500
const int    maxn = 10000;
const int    MAX_I = 100000000;
const double eps = 1e-10;
const double dependence_offset = 2;   //0.5(FDM)     //3.5(ceramic)  //5
const double MAX_D = 1e18;
const double MIN_D = -1e18;
const double dh = 2.0;
const int num_ori_sample = 200; //200

const std::string file_name = ".\\data"; //layer_graph.cpp里面用到了 引用了Mygraph

struct IniData
{
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> sections;
};

inline static std::string Trim(std::string s)
{
	auto notSpace = [](unsigned char ch) { return !std::isspace(ch); };

	s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
	s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
	return s;
}

inline static std::vector<std::string> SplitCsv(const std::string& s)
{
	std::vector<std::string> result;
	std::stringstream ss(s);
	std::string item;

	while (std::getline(ss, item, ',')) {
		item = Trim(item);
		if (!item.empty()) {
			result.push_back(item);
		}
	}

	return result;
}

inline static IniData LoadIni(const std::filesystem::path& iniPath)
{
	IniData data;
	std::ifstream file(iniPath);
	if (!file) {
		std::cerr << "无法打开 ini 文件: " << iniPath << std::endl;
		return data;
	}

	std::string line;
	std::string currentSection;

	while (std::getline(file, line)) {
		line = Trim(line);

		if (line.empty() || line[0] == ';' || line[0] == '#') {
			continue;
		}

		if (line.front() == '[' && line.back() == ']') {
			currentSection = Trim(line.substr(1, line.size() - 2));
			continue;
		}

		const auto pos = line.find('=');
		if (pos == std::string::npos || currentSection.empty()) {
			continue;
		}

		std::string key = Trim(line.substr(0, pos));
		std::string value = Trim(line.substr(pos + 1));

		auto commentPos = value.find_first_of(";#");
		if (commentPos != std::string::npos) {
			value = Trim(value.substr(0, commentPos));
		}

		data.sections[currentSection][key] = value;
	}

	return data;
}

inline static std::string GetIniString(const IniData& ini,
	const std::string& section,
	const std::string& key,
	const std::string& defaultValue = "")
{
	auto secIt = ini.sections.find(section);
	if (secIt == ini.sections.end()) {
		return defaultValue;
	}

	auto keyIt = secIt->second.find(key);
	if (keyIt == secIt->second.end()) {
		return defaultValue;
	}

	return keyIt->second;
}

inline static std::vector<std::string> GetIniStringList(const IniData& ini,
	const std::string& section,
	const std::string& key)
{
	return SplitCsv(GetIniString(ini, section, key));
}

#endif