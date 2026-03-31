#include "Visualization.h"
#include <igl/readOBJ.h>
#include <fstream>

namespace vasco
{
    void createRedBalls(const std::string& file_name,
                        const std::vector<Eigen::MatrixXd>& vis_points)
    {
        Eigen::MatrixXd V_2;
        Eigen::MatrixXi F_2;

        if (!igl::readOBJ("ball.obj", V_2, F_2)) {
			std::cerr << "[createRedBalls] Failed to read ball.obj" << std::endl;
			return;
        }
        std::ofstream all_balls(file_name + "_unaccessivle_pointsxxx.obj");
        for (int i = 0; i < vis_points.size(); i++) {
            for (int j = 0; j < V_2.rows(); j++)
                all_balls << "v "
                << V_2(j, 0) + vis_points[i](0, 0) << " "
                << V_2(j, 1) + vis_points[i](1, 0) << " "
                << V_2(j, 2) + vis_points[i](2, 0) << " "
                << 0.9 << " " << 0.05 << " " << 0.05 << "\n";
            for (int j = 0; j < F_2.rows(); j++)
                all_balls << "f "
                << F_2(j, 0) + i * V_2.rows() + 1 << " "
                << F_2(j, 1) + i * V_2.rows() + 1 << " "
                << F_2(j, 2) + i * V_2.rows() + 1 << "\n";
        }
        all_balls.close();
    }

    void createGreenBalls(const std::string& file_name,
                          const std::vector<Eigen::MatrixXd>& vis_points,
                          const std::vector<double>& color_map)
    {
        Eigen::MatrixXd V_2;
        Eigen::MatrixXi F_2;
        if (!igl::readOBJ("ball.obj", V_2, F_2))
        {
            return;
        }

        std::ofstream all_balls(file_name + "_covering_pointsyyy.obj");
        if (!all_balls.is_open())
        {
            return;
        }

        for (int i = 0; i < static_cast<int>(vis_points.size()); ++i)
        {
            double g = (i < static_cast<int>(color_map.size())) ? color_map[i] * 1.0 : 0.5;
            for (int j = 0; j < V_2.rows(); ++j)
            {
                all_balls << "v "
                          << V_2(j, 0) + vis_points[i](0, 0) << " "
                          << V_2(j, 1) + vis_points[i](1, 0) << " "
                          << V_2(j, 2) + vis_points[i](2, 0) << " "
                          << 0.05 << " " << g << " " << 0.05 << "\n";
            }
            for (int j = 0; j < F_2.rows(); ++j)
            {
                all_balls << "f "
                          << F_2(j, 0) + i * V_2.rows() + 1 << " "
                          << F_2(j, 1) + i * V_2.rows() + 1 << " "
                          << F_2(j, 2) + i * V_2.rows() + 1 << "\n";
            }
        }
        all_balls.close();
    }


	void visualize_layers(const std::string& file_name,
		const Eigen::Vector3f& vn,
		bool judge_continue_additive,
		int id_continue)
	{
		cout << "aa" << endl;
		string input_file_name;
		if (judge_continue_additive == false)
			input_file_name = file_name + ".txt";
		else
			input_file_name = file_name + "_" + to_string(id_continue) + "_subblock.txt";
		cout << input_file_name << endl;
		ifstream ifile(input_file_name);
		int num_contours;
		ifile >> num_contours;
		vector<vector<Eigen::Vector3f>> inputVertices;
		for (int i = 0; i < num_contours; i++) {
			vector<Eigen::Vector3f> currentContour;
			int num_points;
			ifile >> num_points;
			for (int j = 0; j < num_points; j++)
			{
				//cout << "a";
				Eigen::Vector3f currentP;
				ifile >> currentP.x() >> currentP.y() >> currentP.z();
				currentContour.push_back(currentP);
			}
			inputVertices.push_back(currentContour);
		}
		cout << "bb" << num_contours << endl;
		//输入参数
		double R = 0.3;	//切向圆的半径
		int nr = 20;		//用几边形来拟合横截面圆
		int type = 1;	//划分成网格时的类型 0为四边形网格,暂不用 
		//        b---d           b---d 
		// mode 1 | / |    mode 2 | \ |  3为1和2交替 4为2和1交替
		//        a---c           a---c 
		//默认type=1

		double alpha = 1, beta = 1;			//椭圆,alpha*副法向B, beta*正法向N
		Eigen::Vector3f outputColor(1, 1, 255);		//模型颜色,无颜色为(-1,-1,-1)
		string objPath;
		if (judge_continue_additive == false)
			objPath = file_name + "_vis_small.obj";	//输出obj的路径
		else
			objPath = file_name + "_" + to_string(id_continue) + "_subblock_vis_small.obj";	//输出obj的路径


		//调用函数
		vector<Eigen::Vector3f> outputVertices;
		vector<core::Tri3> outputFace;
		//Eigen::Vector3f v1(-10.886800, 8.700820, 64.546097), v2(-13.794400, 9.800380, 62.789299), v3(-22.245300, 2.629520, 68.324501);
		//double na = (v2.y - v1.y) * (v3.z - v1.z) - (v2.z - v1.z) * (v3.y - v1.y);
		//double nb = (v2.z - v1.z) * (v3.x - v1.x) - (v2.x - v1.x) * (v3.z - v1.z);
		//double nc = (v2.x - v1.x) * (v3.y - v1.y) - (v2.y - v1.y) * (v3.x - v1.x);
		//Eigen::Vector3f vn(na, nb, nc);
		//double t = sqrt(vn.x * vn.x + vn.y * vn.y + vn.z * vn.z);
		//vn.x /= t;
		//vn.y /= t;
		//vn.z /= t;

		//Eigen::Vector3f direction(0, 0, 1);	//单位向量
		generateClosedTube(vn, inputVertices, R, nr, type, alpha, beta, outputVertices, outputFace);
		cout << "cc" << endl;
		//输出为obj
		io::outputToObj(outputVertices, outputFace, outputColor, objPath);
		return;
	}

	void visualize_layers_stair_case(const std::string& file_name,
		const std::string& vis_file_contain,
		const Eigen::Vector3f& vn,
		bool judge_continue_additive,
		int id_continue)
	{
		string input_file_name, input_file_name_contain;
		if (judge_continue_additive == false)
			input_file_name = file_name + ".txt";
		else
			input_file_name = file_name + "_" + to_string(id_continue) + "_subblock.txt";


		//输入参数
		double h = 2.0;	//每层的厚度

		double alpha = 1, beta = 1;			//椭圆,alpha*副法向B, beta*正法向N
		Eigen::Vector3f outputColor(1, 1, 255);		//模型颜色,无颜色为(-1,-1,-1)
		string objPath;
		if (judge_continue_additive == false)
			objPath = file_name + "_vis_stair_case.obj";	//输出obj的路径
		else
			objPath = file_name + "_" + to_string(id_continue) + "_subblock_vis_stair_case.obj";	//输出obj的路径

		//调用函数
		vector<Eigen::Vector3f> outputVertices;
		vector<core::Tri3> outputFace;

		//新版读入
		ifstream ifile(input_file_name);
		int num_contours;
		ifile >> num_contours;
		vector<Eigen::Vector3f> contour1, contour2;
		for (int i = 0; i < num_contours; i++) {
			//每个contour,读入是否为环isLoop
			int isLoop = -1;
			ifile >> isLoop;
			contour1.clear();
			contour2.clear();
			//外层contour
			int num_points;
			ifile >> num_points;
			for (int j = 0; j < num_points; j++)
			{
				//cout << "a";
				Eigen::Vector3f currentP;
				ifile >> currentP.x() >> currentP.y() >> currentP.z();
				contour1.push_back(currentP);
			}
			if (isLoop)
			{
				//内层contour
				ifile >> num_points;
				for (int j = 0; j < num_points; j++)
				{
					//cout << "a";
					Eigen::Vector3f currentP;
					ifile >> currentP.x() >> currentP.y() >> currentP.z();
					contour2.push_back(currentP);
				}
			}
			printf("正在处理第%d/%d个轮廓或环.\n", i + 1, num_contours);
			core::newGenerateMesh(vn, isLoop, contour1, contour2, h, outputVertices, outputFace);

		}


		//输出为obj
		io::outputToObj(outputVertices, outputFace, outputColor, objPath);
		return;
	}

}   // namespace vasco


