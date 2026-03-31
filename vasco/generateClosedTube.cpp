#include "generateClosedTube.h"

namespace vasco
{
using vasco::core::Tri3;

bool isZeroClosed(double x)
{
    if (x > -1e-7 && x < 1e-7) {
        return true;
    }
    else {
        return false;
    }
}

bool isEqualClosed(const Eigen::Vector3f& p1, const Eigen::Vector3f& p2)
{
    if (isZeroClosed(p1.x() - p2.x()) && isZeroClosed(p1.y() - p2.y()) && isZeroClosed(p1.z() - p2.z())) {
        return true;
    }
    return false;
}

// double frameAngleClzosed(...) { ... } // 保持注释代码不动

Eigen::Vector3f T1ForOptiClosed, N1ForOptiClosed, T2ForOptiClosed, N2ForOptiClosed, B2ForOptiClosed;

// double functionForOptiClosed(...) { ... } // 保持注释代码不动
// void smoothFramesClosed(...) { ... }      // 保持注释代码不动

void generateFaceClosed(int type,
                        int nr,
                        int numOfPoints,
                        int nstep,
                        int numOfOutputVertices,
                        std::vector<Tri3>& returnFace)
{
    for (int i = 0; i < numOfPoints - 1; i++) {
        std::vector<int> ii;
        for (int j = 0; j < nr; j++) {
            ii.push_back(i * nr + j);
        }
        ii.push_back(i * nr);

        std::vector<int> jj;
        if (i == numOfPoints - 2) {
            int beginIndex = ((nr - nstep) % nr + nr) % nr;
            for (int j = beginIndex; j < nr; j++) {
                jj.push_back(j);
            }
            for (int j = 0; j < nr; j++) {
                jj.push_back(j);
            }
        }
        else {
            for (int j = 0; j < static_cast<int>(ii.size()); j++) {
                jj.push_back(ii[j] + nr);
            }
        }

        for (int k = 0; k < nr; k++) {
            int ia = ii[k], ic = ii[k + 1], ib = jj[k], id = jj[k + 1];
            int currentType = type;
            if (type == 3) {
                currentType = ((i + k) % 2 == 0) ? 1 : 2;
            }
            else if (type == 4) {
                currentType = ((i + k) % 2 == 0) ? 2 : 1;
            }

            if (currentType == 1) {
                returnFace.push_back(Tri3{ ia + numOfOutputVertices, id + numOfOutputVertices, ic + numOfOutputVertices });
                returnFace.push_back(Tri3{ ia + numOfOutputVertices, ib + numOfOutputVertices, id + numOfOutputVertices });
            }
            else if (currentType == 2) {
                returnFace.push_back(Tri3{ ib + numOfOutputVertices, ic + numOfOutputVertices, ia + numOfOutputVertices });
                returnFace.push_back(Tri3{ ib + numOfOutputVertices, id + numOfOutputVertices, ic + numOfOutputVertices });
            }
            else {
                std::cout << "Unexpected type: " << currentType << "!" << std::endl;
            }
        }
    }
}

void generateClosedTube(Eigen::Vector3f direction,
                        std::vector<std::vector<Eigen::Vector3f> > inputVertices,
                        double R,
                        int nr,
                        int type,
                        double alpha,
                        double beta,
                        std::vector<Eigen::Vector3f>& outputVertices,
                        std::vector<Tri3>& outputFace)
{
    int numOfContours = static_cast<int>(inputVertices.size());
    for (int contourIndex = 0; contourIndex < numOfContours; contourIndex++) {
        // 取出当前轮廓顶点序列
        std::vector<Eigen::Vector3f> currentContour = inputVertices[contourIndex];
        int numOfPoints = static_cast<int>(currentContour.size());

        // 20220122 自己计算每个顶点的T,N,B
        // T 切向量 顶点所在两个线段切向之和
        // N 正法向 首尾端点任选一个垂直T的即可 其余顶点为两个线段所成夹角一半
        // B 副法向 B=T×N
        std::vector<Eigen::Vector3f> T;
        std::vector<Eigen::Vector3f> N;
        std::vector<Eigen::Vector3f> B;

        std::vector<Eigen::Vector3f> dp; // 两个顶点的差值
        for (int i = 0; i < static_cast<int>(currentContour.size()) - 1; i++) {
            dp.push_back(currentContour[i + 1] - currentContour[i]);
        }
        for (int i = 0; i < static_cast<int>(dp.size()); i++) {
            float len = dp[i].norm();
            dp[i] = dp[i] / len;
        }
        T.push_back(dp[0]);
        for (int i = 0; i < static_cast<int>(dp.size()) - 1; i++) {
            T.push_back(dp[i] + dp[i + 1]);
        }
        T.push_back(dp[dp.size() - 1]);

        for (int i = 0; i < static_cast<int>(T.size()); i++) {
            float len = T[i].norm();
            T[i] = T[i] / len;
        }

        Eigen::Vector3f TT = T[0] + T[T.size() - 1];
        TT = TT / TT.norm();
        T[0] = TT;
        T[T.size() - 1] = TT;

        int nstep = 0;

        N.resize(T.size());
        B.resize(T.size());
        for (int i = 0; i < static_cast<int>(N.size()); i++) {
            N[i] = direction;
            B[i] = T[i].cross(N[i]);
        }

        int numOfOutputVertices = static_cast<int>(outputVertices.size()); // 加入之前已经有这么多顶点了

        std::vector<float> ss;
        std::vector<float> rr;
        for (int i = 0; i < nr + 1; i++) {
            float tt = 2 * static_cast<float>(M_PI) * i / nr;
            ss.push_back(R * std::cos(tt + static_cast<float>(M_PI) / 4) * static_cast<float>(alpha));
            rr.push_back(R * std::sin(tt + static_cast<float>(M_PI) / 4) * static_cast<float>(beta));
        }
        for (int i = 0; i < numOfPoints - 1; i++) {
            Eigen::Vector3f thePoint = currentContour[i];
            for (int k = 0; k < nr; k++) {
                Eigen::Vector3f p = B[i] * ss[k] + N[i] * rr[k];
                p = thePoint + p;
                outputVertices.push_back(p);
            }
        }

        // 计算面片的顶点连接情况
        generateFaceClosed(type, nr, numOfPoints, nstep, numOfOutputVertices, outputFace);
    }
}
} // namespace vasco