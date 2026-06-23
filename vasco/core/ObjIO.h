#pragma once
#include <vector>
#include <string>
#include <fstream>
#include "GeometryUtils.h"
#include "Types.h"


namespace vasco::io {

	// (-1,-1,-1) 表示不写颜色
	// 返回 true 成功, false 打开文件失败
	inline bool outputToObj(const std::vector<Eigen::Vector3f>& vertices,
		const std::vector<vasco::core::Tri3>& faces,
		const Eigen::Vector3f& color,
		const std::string& filename)
	{
		//检查是否需要输出颜色
		std::ofstream fout(filename, std::ios::out);
		if (!fout.is_open())
			return false;

		const bool writeColor = !(color == Eigen::Vector3f(-1.0, -1.0, -1.0));

		for (const auto& v : vertices) {
			if (writeColor) {
				fout << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << ' '
					<< color.x() << ' ' << color.y() << ' ' << color.z() << '\n';
			}
			else {
				fout << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << '\n';
			}
		}
		for (const auto& f : faces) {
			//if (f.v1 < 0 || f.v2 < 0 || f.v3 < 0) continue;
			fout << "f " << (f[0] + 1) << ' ' << (f[1] + 1) << ' ' << (f[2] + 1) << '\n';
		}
		return true;
	}

	// (-1,-1,-1) 表示不写颜色
	// 返回 true 成功, false 打开文件失败
	inline bool outputContourToObj(const std::vector<std::vector<Eigen::Vector3d>>& contours,
		const Eigen::Vector3d& color,
		const std::string& filename)
	{
		//检查是否需要输出颜色
		std::ofstream fout(filename, std::ios::out);
		if (!fout.is_open())
			return false;

		const bool writeColor = !(color == Eigen::Vector3d(-1.0, -1.0, -1.0));

		int vertexCount = 0;
		for (const auto& c : contours) {
			for (const auto& v : c) {
				if (writeColor) {
					fout << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << ' '
						<< color.x() << ' ' << color.y() << ' ' << color.z() << '\n';
				}
				else {
					fout << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << '\n';
				}
			}
			fout << "f";
			for (int i = 0; i < c.size(); ++i) {
				fout << ' ' << (vertexCount + i + 1);
			}
			fout << '\n';
			vertexCount += c.size();
		}
		return true;
	}

} // namespace vasco::io