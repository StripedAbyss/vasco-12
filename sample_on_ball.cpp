#include "sample_on_ball.h"


void SAMPLE_ON_BALL::OrientationSamplePoints()
{
    //�°��˹����������ƽ�Ƕȣ�
    sample_points.clear();
    sample_points.push_back(Eigen::Vector3d(0, 0, 1));

    double goldenAngle = 3.1415926 * (3.0 - sqrt(5.0));
    double offset = 2.0 / num_ori_sample;

    for (int i = 0; i < num_ori_sample; ++i) {
        double y = i * offset - 1.0 + (offset / 2.0);
        double radiusXY = sqrt(1.0 - y * y);
        double theta = goldenAngle * i;

        double x = cos(theta) * radiusXY;
        double z = sin(theta) * radiusXY;

        if (z > -0.02)
            sample_points.push_back({ x * 1, y * 1, z * 1 });
    }
}

void SAMPLE_ON_BALL::OrientationSamplePoints_2()
{
    //�°��˹����������ƽ�Ƕȣ�
    sample_points.clear();
    sample_points.push_back(Eigen::Vector3d(0, 0, 1));

    double goldenAngle = 3.1415926 * (3.0 - sqrt(5.0));
    double offset = 2.0 / num_ori_sample;

    for (int i = 0; i < num_ori_sample; ++i) {
        double y = i * offset - 1.0 + (offset / 2.0);
        double radiusXY = sqrt(1.0 - y * y);
        double theta = goldenAngle * i;

        double x = cos(theta) * radiusXY;
        double z = sin(theta)    * radiusXY;

        if (z > 0.5)
            sample_points.push_back({ x * 1, y * 1, z * 1 });
    }
}
