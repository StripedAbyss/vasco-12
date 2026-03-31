#pragma once
#include "Types.h"
#include "CgalUtils.h"
#include <opencv2/core.hpp>
#include <cmath>
#include "../../PolygonIntersection.h"
#include "../../layer_graph.h"
#include <Eigen/Dense>

namespace vasco::core {

	inline double distance3d(const cv::Point3d& a, const cv::Point3d& b) noexcept {
		double distance = sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
		return distance;
	}
	inline double distance3d(const Vec3& a, const Vec3& b) noexcept {
		double distance = sqrt((a[0] - b[0]) * (a[0] - b[0]) +
			(a[1] - b[1]) * (a[1] - b[1]) +
			(a[2] - b[2]) * (a[2] - b[2]));
		return distance;
	}
	inline double distance2d(const cv::Point2d& a, const cv::Point2d& b) noexcept {
		double distance = sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
		return distance;
	}
	inline double distanceVec3(const Vec3& a, const Vec3& b) noexcept {
		double distance = sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2));
		return distance;
	}
	/*inline Vec3 faceNormal(const Vec3& v1, const Vec3& v2, const Vec3& v3) noexcept {
		return {
			(v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]),
			(v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]),
			(v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0])
		};
	}*/
	inline Eigen::Vector3d faceNormal(const Vec3& v1, const Vec3& v2, const Vec3& v3) noexcept {
		//double ans = (v2[0] - v1[0]) * (v2[1] - v3[1]) - (v2[1] - v1[1]) * (v2[0] - v3[0]);
		//if (ans > 0)	//is clockwise
		//	swap(v2, v3);

		double na = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]);
		double nb = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]);
		double nc = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]);

		Eigen::Vector3d vn(na, nb, nc);
		return vn;
	}

	//判断两个向量几乎相同
	inline bool isEqual(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2)
	{
		// 使用 Eigen 的 L-infinity 范数代替，逻辑与原代码完全一致（检查最大分量差值）
		// (v1 - v2).lpNorm<Eigen::Infinity>() 等价于 max(|dx|, |dy|, |dz|)
		return (v1 - v2).lpNorm<Eigen::Infinity>() < 1e-2f;

		/*if (v1.x() - v2.x() > -1e-2 && v1.x() - v2.x() < 1e-2)
			if (v1.y() - v2.y() > -1e-2 && v1.y() - v2.y() < 1e-2)
				if (v1.z() - v2.z() > -1e-2 && v1.z() - v2.z() < 1e-2)
					return true;
		return false;*/
	}

	//判断两个向量方向相同
	//因为经常会有法向跟打印方向差了0.0几导致不相等,所以写个这个
	inline bool isSameDirection(Eigen::Vector3f v1, Eigen::Vector3f v2)
	{
		if (v1.dot(v2) > 0)
			return true;
		else
			return false;
	}

	inline double triangleArea(const Vec3& v1, const Vec3& v2, const Vec3& v3) noexcept {

		Eigen::Vector3d vec1, vec2;
		vec1.x() = v2[0] - v1[0];
		vec1.y() = v2[1] - v1[1];
		vec1.z() = v2[2] - v1[2];
		vec2.x() = v3[0] - v1[0];
		vec2.y() = v3[1] - v1[1];
		vec2.z() = v3[2] - v1[2];

		double area = vec1.cross(vec2).norm() / 2;
		return area;
	}

	//计算三角形有向面积
	inline double calculateTriArea(const Eigen::Vector3f& direction, const Eigen::Vector3f& p1, const Eigen::Vector3f& p2, const Eigen::Vector3f& p3)
	{
		// 计算三角形法向量 (叉乘)
		Eigen::Vector3f N = (p2 - p1).cross(p3 - p1);
		
		// 计算面积 (模长的一半)
		double S = 0.5 * N.norm();

		// 判断方向 (点乘 > 0 表示方向相同)
		if (N.dot(direction) > 0)
			return S;
		else
			return -S;
	}

	inline int isAnticlockwise(cv::Point3f p1, cv::Point3f p2, cv::Point3f p3)
	{
		//
		double ans = (p2.x - p1.x) * (p2.y - p3.y) - (p2.y - p1.y) * (p2.x - p3.x);
		//std::cout << ans << '\n';
		if (ans < 0)		//is anticlockwise
			return 1;
		else if (ans == 0)	//is line
			return 2;
		else if (ans > 0)			//is clockwise
			return 3;
		return -1;
	}

	bool checkHaveNoOtherPoint(const cv::Point3f& p1,
		const cv::Point3f& p2,
		const cv::Point3f& p3,
		const std::list<cv::Point3f>& pointList);

	bool checkHaveNoOtherPoint(const Eigen::Vector3f& direction,
		const Eigen::Vector3f& p1,
		const Eigen::Vector3f& p2,
		const Eigen::Vector3f& p3,
		const std::list<Eigen::Vector3f>& pointList);

	int isAnticlockwise(const Eigen::Vector3f& direction, const Eigen::Vector3f& p1, const Eigen::Vector3f& p2, const Eigen::Vector3f& p3);
	bool pointsAreAnticlockwise(Eigen::Vector3f direction, std::vector<Eigen::Vector3f> inputVertices);

	void polygonToTriangle(const Eigen::Vector3f& direction,
		const std::vector<Eigen::Vector3f>& inputVertices,
		bool isAnti,
		double h,
		std::vector<Eigen::Vector3f>& outputVertices,
		std::vector<Tri3>& outputFace);

	void generateWall(Eigen::Vector3f direction, std::vector<Eigen::Vector3f> inputVertices, bool isAnti, double h, std::vector<Eigen::Vector3f>& outputVertices, std::vector<Tri3>& outputFace);

	void generateLoop(Eigen::Vector3f direction, std::vector<Eigen::Vector3f> contour1, std::vector<Eigen::Vector3f> contour2, double h, std::vector<Eigen::Vector3f>& outputVertices, std::vector<Tri3>& outputFace);

	void newGenerateMesh(Eigen::Vector3f direction, bool isLoop, std::vector<Eigen::Vector3f> contour1, std::vector<Eigen::Vector3f> contour2, double h, std::vector<Eigen::Vector3f>& outputVertices, std::vector<Tri3>& outputFace);

	void generateMesh(Eigen::Vector3f direction, std::vector<std::vector<Eigen::Vector3f>> inputVertices, double h, std::vector<Eigen::Vector3f>& outputVertices, std::vector<Tri3>& outputFace);

	void getMeshNormal(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, Eigen::MatrixXd& normals);

	bool isPointInTriangle(const Eigen::Vector3d& pt, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3);

	void vec3toMatrix(const std::vector<Vec3>& inputVerts, Eigen::MatrixXd& vertices);

	void tri3ToMatrix(const std::vector<Tri3>& inputTris, Eigen::MatrixXi& faces);

} // namespace vasco::core



//Vector3 HybridManufacturing::calculate_normal(VEctor v1, VEctor v2, VEctor v3)
//{
//	//double ans = (v2[0] - v1[0]) * (v2[1] - v3[1]) - (v2[1] - v1[1]) * (v2[0] - v3[0]);
//	//if (ans > 0)	//is clockwise
//	//	swap(v2, v3);
//
//	double na = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]);
//	double nb = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]);
//	double nc = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]);
//
//	Vector3 vn(na, nb, nc);
//	return vn;
//}
//
//double HybridManufacturing::calculate_area(VEctor v1, VEctor v2, VEctor v3)
//{
//	Eigen::Vector3d vec1, vec2;
//	vec1.x() = v2[0] - v1[0];
//	vec1.y() = v2[1] - v1[1];
//	vec1.z() = v2[2] - v1[2];
//	vec2.x() = v3[0] - v1[0];
//	vec2.y() = v3[1] - v1[1];
//	vec2.z() = v3[2] - v1[2];
//
//	double area = vec1.cross(vec2).norm() / 2;
//	return area;
//}