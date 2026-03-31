#include <cmath>
#include <iostream>
#include <vector>
#include "core/GeometryUtils.h"
#include "core/Types.h"

#define M_PI 3.14159265358979323846

namespace vasco
{
// 计算两个帧之间的旋转角度
double frameAngleClosed(Eigen::Vector3f T1, Eigen::Vector3f N1, Eigen::Vector3f T2, Eigen::Vector3f N2);

// 对帧进行平滑
void smoothFramesClosed(std::vector<Eigen::Vector3f> T,
                        std::vector<Eigen::Vector3f>& N,
                        std::vector<Eigen::Vector3f>& B,
                        int nr,
                        int& nstep);

// 计算面片,其实这里用不到nstep
void generateFaceClosed(int type,
                        int nr,
                        int numOfPoints,
                        int nstep,
                        int numOfOutputVertices,
                        std::vector<vasco::core::Tri3>& outputFace);

// 输入/输出接口
void generateClosedTube(Eigen::Vector3f direction,
                        std::vector<std::vector<Eigen::Vector3f> > inputVertices,
                        double R,
                        int nr,
                        int type,
                        double alpha,
                        double beta,
                        std::vector<Eigen::Vector3f>& outputVertices,
                        std::vector<vasco::core::Tri3>& outputFace);
} // namespace vasco

//输入
//vector<vector<Point3ff> > inputVertices;//二维数组 i条轮廓 每条闭合轮廓<a,b,...,g,a>
//double R = 1.0;				//切向圆的半径
//int nr = 20;					//用几边形来拟合切向圆
//int type = 1;					//划分成网格时的类型 0为四边形网格,暂不用 
//        b---d           b---d 
// mode 1 | / |    mode 2 | \ |  3为1和2交替 4为2和1交替
//        a---c           a---c 
//double alpha = 1, beta = 1;	//alpha*副法向B, beta*正法向N

//输出
//vector<Point3ff> outputVertices;	//模型顶点序列
//vector<face> outputFace;			//面片的顶点索引

//trivial 1:输入的采样点不能相邻三点在同一条线上,需要增加输入检查
//trivial 2:求解RMF时dlib库可能有bug,写了个cout
