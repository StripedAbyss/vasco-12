#include "GeometryUtils.h"

namespace vasco::core {

	bool checkHaveNoOtherPoint(const cv::Point3f &p1, 
		const cv::Point3f &p2, 
		const cv::Point3f &p3, 
		const std::list<cv::Point3f> &pointList)
	{
		for (const auto& p : pointList)
		{
			if (p == p1 || p == p2 || p == p3)
				continue;
			// if point p is in the triangle, (p,p1,p2),(p,p2,p3),(p,p3,p1) should be anticlockwise
			if (isAnticlockwise(p, p1, p2) == 1 &&
				isAnticlockwise(p, p2, p3) == 1 &&
				isAnticlockwise(p, p3, p1) == 1)
				return false;
		}
		return true;
	}

	//检测三点中间是否还有其他顶点
	bool checkHaveNoOtherPoint(const Eigen::Vector3f &direction, 
		const Eigen::Vector3f &p1, 
		const Eigen::Vector3f &p2, 
		const Eigen::Vector3f &p3, 
		const std::list<Eigen::Vector3f> &pointList)
	{
		for (const auto& p : pointList) {
			if (p == p1 || p == p2 || p == p3) {
				continue;
			}
			if (isAnticlockwise(direction, p, p1, p2) == 1 &&
				isAnticlockwise(direction, p, p2, p3) == 1 &&
				isAnticlockwise(direction, p, p3, p1) == 1) {
				return false;
			}
		}
		return true;
	}

	//在-direction方向下看去,三点是否为逆时针
	int isAnticlockwise(const Eigen::Vector3f& direction,
		const Eigen::Vector3f& p1,
		const Eigen::Vector3f& p2,
		const Eigen::Vector3f& p3)
	{
		// p2p1 叉乘 p2p3 与 direction 方向相同 则为逆时针
		// 注意：原代码逻辑是 "与 direction * (-1) 相同则为逆时针"，即与 -direction 同向

		Eigen::Vector3f p21 = (p1 - p2).normalized();
		Eigen::Vector3f p23 = (p3 - p2).normalized();

		// 计算叉积
		Eigen::Vector3f ans = p21.cross(p23);
		double ansNorm2 = ans.squaredNorm();

		// 判断是否共线 (sin(theta) 非常小)
		// sin(1度) ≈ 0.017，原代码阈值为 0.02
		if (ansNorm2 < 0.02 * 0.02)
			return 2; // is line

		// 判断方向
		// 原代码: isSameDirection(ans, direction * (-1)) -> ans . (-direction) > 0 -> ans . direction < 0
		if (ans.dot(direction) < 0)
			return 1; // is anticlockwise
		else
			return 3; // is clockwise
	}

	//判断这一层顶点在-direction方向下看去是否为逆时针
	bool pointsAreAnticlockwise(Eigen::Vector3f direction, std::vector<Eigen::Vector3f> inputVertices)
	{
		//输入检查
		const std::size_t n = inputVertices.size();
		if (n < 3) {
			return false;
		}
		//计算有向面积
		double S = 0;
		const Eigen::Vector3f &p1 = inputVertices[0];
		for (std::size_t i = 1; i + 1 < n; ++i) 
		{
			const Eigen::Vector3f &p2 = inputVertices[i];
			const Eigen::Vector3f &p3 = inputVertices[i + 1];
			S = S + calculateTriArea(direction, p1, p2, p3);
		}
		if (S > 0)
			return true;
		else
			return false;

	}


	//任意(凹)多边形转换为多个三角形
	//顶面在-direction方向下看去是逆时针
	//底面在direction方向下看去是逆时针
	//inputVertices此时为一层多边形轮廓
	void polygonToTriangle(const Eigen::Vector3f& direction,
		const std::vector<Eigen::Vector3f>& inputVertices,
		bool isAnti,
		double h,
		std::vector<Eigen::Vector3f>& outputVertices,
		std::vector<Tri3>& outputFace)
	{
		const int numOfPoint = static_cast<int>(inputVertices.size()) - 1; // 最后一个与首点相同
		if (numOfPoint < 3) {
			return;
		}

		// 按需要的朝向、并上移 h
		std::vector<Eigen::Vector3f> verts;
		verts.reserve(numOfPoint);
		if (!isAnti) {
			for (int i = numOfPoint - 1; i >= 0; --i) {
				verts.push_back(inputVertices[i] + direction * static_cast<float>(h));
			}
		} else {
			for (int i = 0; i < numOfPoint; ++i) {
				verts.push_back(inputVertices[i] + direction * static_cast<float>(h));
			}
		}

		// 索引循环队列
		std::vector<int> idx(numOfPoint);
		std::iota(idx.begin(), idx.end(), 0);

		auto isEar = [&](int a, int b, int c) -> bool {
			// 必须是逆时针耳
			if (isAnticlockwise(direction, verts[a], verts[b], verts[c]) != 1) {
				return false;
			}
			// 检查是否有其他点在三角形内
			for (int v : idx) {
				if (v == a || v == b || v == c) {
					continue;
				}
				if (isAnticlockwise(direction, verts[v], verts[a], verts[b]) == 1 &&
					isAnticlockwise(direction, verts[v], verts[b], verts[c]) == 1 &&
					isAnticlockwise(direction, verts[v], verts[c], verts[a]) == 1) {
					return false;
				}
			}
			return true;
		};

		while (static_cast<int>(idx.size()) >= 3) {
			bool clipped = false;
			const int m = static_cast<int>(idx.size());
			for (int k = 0; k < m; ++k) {
				int ia = idx[k];
				int ib = idx[(k + 1) % m];
				int ic = idx[(k + 2) % m];

				if (!isEar(ia, ib, ic)) {
					continue;
				}

				const int base = static_cast<int>(outputVertices.size());
				outputVertices.push_back(verts[ia]);
				outputVertices.push_back(verts[ib]);
				outputVertices.push_back(verts[ic]);
				outputFace.push_back(Tri3{ base, base + 1, base + 2 });

				// 剪掉耳朵
				idx.erase(idx.begin() + ((k + 1) % m));
				clipped = true;
				break;
			}
			if (!clipped) {
				// 退化或数值问题，防止死循环
				break;
			}
		}
	}


	//生成侧壁
	void generateWall(Eigen::Vector3f direction, std::vector<Eigen::Vector3f> inputVertices, bool isAnti, double h, std::vector<Eigen::Vector3f>& outputVertices, std::vector<Tri3>& outputFace)
	{
		//对每一对顶点x和y,生成打印方向对应的xx和yy,然后连接成两个三角形面片
		for (int xIndex = 0; xIndex < inputVertices.size() - 1; xIndex++)
		{
			//构造顶点
			Eigen::Vector3f x = inputVertices[xIndex];
			Eigen::Vector3f y = inputVertices[xIndex + 1];
			Eigen::Vector3f xx = x + direction * h;
			Eigen::Vector3f yy = y + direction * h;
			int indexOfX = outputVertices.size();
			outputVertices.push_back(x);		//indexOfX
			outputVertices.push_back(xx);	//indexOfX+1
			outputVertices.push_back(y);		//indexOfX+2
			outputVertices.push_back(yy);	//indexOfX+3

			if (isAnti)
			{
				//构造面片[x,y,xx]和[y,yy,xx]
				outputFace.push_back(Tri3{ indexOfX, indexOfX + 2, indexOfX + 1 });
				outputFace.push_back(Tri3{ indexOfX + 2, indexOfX + 3, indexOfX + 1 });
			}
			else
			{
				//构造面片[x,xx,y]和[y,xx,yy]
				outputFace.push_back(Tri3{ indexOfX, indexOfX + 1, indexOfX + 2 });
				outputFace.push_back(Tri3{ indexOfX + 2, indexOfX + 1, indexOfX + 3 });
			}
		}

	}

	//生成环
	//顶面在-direction方向下看去是逆时针
	//底面在direction方向下看去是逆时针
	void generateLoop(Eigen::Vector3f direction, std::vector<Eigen::Vector3f> contour1, std::vector<Eigen::Vector3f> contour2, double h, std::vector<Eigen::Vector3f>& outputVertices, std::vector<Tri3>& outputFace)
	{
		//这个函数里会默认数组首尾顶点是不同的,所以需要删掉最后一个顶点
		contour1.pop_back();
		contour2.pop_back();
		//确认起始点,按最近策略连接顶点
		double dis = (contour1[0]-contour2[0]).squaredNorm();
		int indexOfVertex = 0;	//挑一个离contour1[0]最近的顶点
		for (int i = 1; i < contour2.size(); i++)
		{
			double dis_i = (contour1[0] - contour2[i]).squaredNorm();
			if (dis > dis_i)
			{
				dis = dis_i;
				indexOfVertex = i;
			}
		}
		//重新排序contour2,此时两个起始点最近
		std::vector<Eigen::Vector3f> temp;
		for (int i = 0; i < contour2.size(); i++)
			temp.push_back(contour2[(i + indexOfVertex) % contour2.size()]);
		contour2 = temp;

		//循环:每次拿出一个contour1和一个contour2的顶点,并在contour1[i+1]和contour2[j+1]这两个候选顶点中选择一个,连接为三角形
		//选择方法:数目差/组成的三角形内部不包括候选顶点和环中心/三点不在一条直线上
		//退出条件:有个三角形的顶点是初始两个顶点之一,即连接了一圈回来了
		int index1 = 0, index2 = 0;
		while (index1 != contour1.size() && index2 != contour2.size())
		{
			Eigen::Vector3f p1 = contour1[index1];
			Eigen::Vector3f p2 = contour2[index2];
			Eigen::Vector3f p3 = contour1[(index1 + 1) % contour1.size()];
			bool chooseContour1 = true;	//p3选contour1上的
			//首先比较剩余顶点的数目,优先选剩的多的
			int numOfV1 = contour1.size() - index1, numOfV2 = contour2.size() - index2;
			if (numOfV1 < numOfV2)
			{
				p3 = contour2[(index2 + 1) % contour2.size()];
				chooseContour1 = false;
			}
			//接下来判断该三角形内部不包含候选顶点和环中心
			bool check = true;

			if (chooseContour1)	//外-内-外应该是顺时针,且不包含contour2[index2+1],且该点和p1p3,p2p3不在一条直线上
			{
				if (isAnticlockwise(direction, p1, p2, p3) != 3)
					check = false;
				Eigen::Vector3f p = contour2[(index2 + 1) % contour2.size()];
				if (isAnticlockwise(direction, p, p2, p1) == 1 && isAnticlockwise(direction, p, p3, p2) == 1 && isAnticlockwise(direction, p, p1, p3) == 1)
					check = false;
				if (isAnticlockwise(direction, p1, p, p3) == 2 || isAnticlockwise(direction, p2, p, p3) == 2)
					check = false;
			}
			else //外-内-内应该是顺时针,且不包含contour1[index1+1],且该点和p1p3,p2p3不在一条直线上
			{
				if (isAnticlockwise(direction, p1, p2, p3) != 3)
					check = false;
				Eigen::Vector3f p = contour1[(index1 + 1) % contour1.size()];
				if (isAnticlockwise(direction, p, p2, p1) == 1 && isAnticlockwise(direction, p, p3, p2) == 1 && isAnticlockwise(direction, p, p1, p3) == 1)
					check = false;
				if (isAnticlockwise(direction, p1, p, p3) == 2 || isAnticlockwise(direction, p2, p, p3) == 2)
					check = false;
			}

			//总能选出一个来吧,就不加两个都不能选的特殊处理了
			if (check == false && chooseContour1 == true)
			{
				chooseContour1 = false;
				p3 = contour2[(index2 + 1) % contour2.size()];
			}
			else if (check == false && chooseContour1 == false)
			{
				chooseContour1 = true;
				p3 = contour1[(index1 + 1) % contour1.size()];
			}

			//顶面:将该面片[p1,p3,p2]+direction*h (逆时针)添加到数组中
			int indexOfP = outputVertices.size();
			outputVertices.push_back(p1 + direction * h);
			outputVertices.push_back(p3 + direction * h);
			outputVertices.push_back(p2 + direction * h);
			outputFace.push_back(Tri3{ indexOfP, indexOfP + 1, indexOfP + 2 });

			//底面:将该面片[p1,p2,p3] (顺时针)添加到数组中
			indexOfP = outputVertices.size();
			outputVertices.push_back(p1);
			outputVertices.push_back(p2);
			outputVertices.push_back(p3);
			outputFace.push_back(Tri3{ indexOfP, indexOfP + 1, indexOfP + 2 });

			//这里不对数组大小取余
			if (chooseContour1)
				index1++;
			else
				index2++;


		}
		//最后的特殊处理,此时有一侧contour连接完毕,将另一侧连接过来
		if (index1 == contour1.size())
		{
			while (index2 != contour2.size())
			{
				Eigen::Vector3f p1 = contour1[0];
				Eigen::Vector3f p2 = contour2[index2];
				Eigen::Vector3f p3 = contour2[(index2 + 1) % contour2.size()];

				//应该不需要判断能否加入吧

				//顶面:将该面片[p1,p3,p2]+direction*h (逆时针)添加到数组中
				int indexOfP = outputVertices.size();
				outputVertices.push_back(p1 + direction * h);
				outputVertices.push_back(p3 + direction * h);
				outputVertices.push_back(p2 + direction * h);
				outputFace.push_back(Tri3{ indexOfP, indexOfP + 1, indexOfP + 2 });

				//底面:将该面片[p1,p2,p3] (顺时针)添加到数组中
				indexOfP = outputVertices.size();
				outputVertices.push_back(p1);
				outputVertices.push_back(p2);
				outputVertices.push_back(p3);
				outputFace.push_back(Tri3{ indexOfP, indexOfP + 1, indexOfP + 2 });

				index2++;
			}
		}
		else if (index2 == contour2.size())
		{
			while (index1 != contour1.size())
			{
				Eigen::Vector3f p1 = contour1[index1];
				Eigen::Vector3f p2 = contour2[0];
				Eigen::Vector3f p3 = contour1[(index1 + 1) % contour1.size()];

				//应该不需要判断能否加入吧

				//顶面:将该面片[p1,p3,p2]+direction*h (逆时针)添加到数组中
				int indexOfP = outputVertices.size();
				outputVertices.push_back(p1 + direction * h);
				outputVertices.push_back(p3 + direction * h);
				outputVertices.push_back(p2 + direction * h);
				outputFace.push_back(Tri3{ indexOfP, indexOfP + 1, indexOfP + 2 });

				//底面:将该面片[p1,p2,p3] (顺时针)添加到数组中
				indexOfP = outputVertices.size();
				outputVertices.push_back(p1);
				outputVertices.push_back(p2);
				outputVertices.push_back(p3);
				outputFace.push_back(Tri3{ indexOfP, indexOfP + 1, indexOfP + 2 });

				index1++;
			}
		}

	}


	//新-生成面片 支持环结构
	void newGenerateMesh(Eigen::Vector3f direction, bool isLoop, std::vector<Eigen::Vector3f> contour1, std::vector<Eigen::Vector3f> contour2, double h, std::vector<Eigen::Vector3f>& outputVertices, std::vector<Tri3>& outputFace)
	{
		//isLoop为false,和原来一致,生成侧壁和顶面底面的三角面片
		//此时contour2不使用
		if (!isLoop)
		{
			//判断顶点方向
			bool isAnti = pointsAreAnticlockwise(direction, contour1);
			//生成侧壁
			generateWall(direction, contour1, isAnti, h, outputVertices, outputFace);
			//顶面传入h
			polygonToTriangle(direction, contour1, isAnti, h, outputVertices, outputFace);
			//底面传入h=0
			polygonToTriangle(direction * (-1), contour1, !isAnti, 0, outputVertices, outputFace);
		}
		else //isLoop为true,为环结构,contour1为外层
		{
			//有可能两层环方向不同,将其都变为在-direction方向下看去是逆时针
			if (!pointsAreAnticlockwise(direction, contour1))
			{
				std::vector<Eigen::Vector3f> temp;
				for (int i = contour1.size() - 1; i >= 0; i--)
					temp.push_back(contour1[i]);
				contour1 = temp;
			}
			if (!pointsAreAnticlockwise(direction, contour2))
			{
				std::vector<Eigen::Vector3f> temp;
				for (int i = contour2.size() - 1; i >= 0; i--)
					temp.push_back(contour2[i]);
				contour2 = temp;
			}
			bool isAnti = true;
			//生成侧壁
			generateWall(direction, contour1, isAnti, h, outputVertices, outputFace);
			generateWall(direction, contour2, !isAnti, h, outputVertices, outputFace);	//注意法向是反的,所以输入了!isAnti
			//生成顶面和地面
			generateLoop(direction, contour1, contour2, h, outputVertices, outputFace);
		}
	}


	//生成侧壁和顶面底面的三角面片
	void generateMesh(Eigen::Vector3f direction, std::vector<std::vector<Eigen::Vector3f> > inputVertices, double h, std::vector<Eigen::Vector3f>& outputVertices, std::vector<Tri3>& outputFace)
	{
		//对于layerIndex层,生成侧壁,底面和顶面
		for (int layerIndex = 0; layerIndex < inputVertices.size(); layerIndex++)
		{
			//判断顶点方向
			bool isAnti = pointsAreAnticlockwise(direction, inputVertices[layerIndex]);
			//生成侧壁
			generateWall(direction, inputVertices[layerIndex], isAnti, h, outputVertices, outputFace);
			//顶面传入h
			polygonToTriangle(direction, inputVertices[layerIndex], isAnti, h, outputVertices, outputFace);
			//底面传入h=0
			polygonToTriangle(direction * (-1), inputVertices[layerIndex], !isAnti, 0, outputVertices, outputFace);
			printf("正在处理%d/%d.\n", layerIndex + 1, (int)inputVertices.size());
		}
	}

	void getMeshNormal(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, Eigen::MatrixXd& normals)
	{
		normals.resize(faces.rows(), 3);
		for (int i = 0; i < faces.rows(); i++)
		{
			Eigen::Vector3d v1, v2, v3;
			v1 = vertices.row(faces(i, 0));
			v2 = vertices.row(faces(i, 1));
			v3 = vertices.row(faces(i, 2));
			Eigen::Vector3d vec1, vec2;
			vec1 = v2 - v1;
			vec2 = v3 - v1;
			Eigen::Vector3d normal = vec1.cross(vec2).normalized();
			normals.row(i) = normal;
		}
	}

	bool isPointInTriangle(const Eigen::Vector3d& pt, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3) {
		const double eps_here = 0.003;
		Eigen::Vector3d v2v1 = v2 - v1;
		Eigen::Vector3d v3v1 = v3 - v1;
		Eigen::Vector3d ptv1 = pt - v1;
		double d00 = v2v1.dot(v2v1);
		double d01 = v2v1.dot(v3v1);
		double d11 = v3v1.dot(v3v1);
		double d20 = ptv1.dot(v2v1);
		double d21 = ptv1.dot(v3v1);
		double denom = d00 * d11 - d01 * d01;
		if (denom == 0) {
			return false; // Degenerate triangle
		}
		double v = (d11 * d20 - d01 * d21) / denom;
		double w = (d00 * d21 - d01 * d20) / denom;
		double u = 1.0 - v - w;
		return (u >= -eps_here) && (v >= -eps_here) && (w >= -eps_here);
	}

	void vec3toMatrix(const std::vector<Vec3>& inputVerts, Eigen::MatrixXd& vertices)
	{
		vertices.resize(inputVerts.size(), 3);
		for (int i = 0; i < inputVerts.size(); i++)
		{
			vertices(i, 0) = inputVerts[i][0];
			vertices(i, 1) = inputVerts[i][1];
			vertices(i, 2) = inputVerts[i][2];
		}
	}

	void tri3ToMatrix(const std::vector<Tri3>& inputTris, Eigen::MatrixXi& faces)
	{
		faces.resize(inputTris.size(), 3);
		for (int i = 0; i < inputTris.size(); i++)
		{
			faces(i, 0) = inputTris[i][0];
			faces(i, 1) = inputTris[i][1];
			faces(i, 2) = inputTris[i][2];
		}
	}

} // namespace vasco::core