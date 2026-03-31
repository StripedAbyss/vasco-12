#pragma once
#include<cstdio>
#include<iostream>
#include<algorithm>
#include<cstring>
#include<cmath>
using namespace std;
#define maxn_poly 510

//inline int sig(double d) {
//    return(d > eps) - (d < -eps);
//}
struct MyPoint {
    double x, y; 
    MyPoint() { x = 0; y = 0; }
    MyPoint(double x, double y) :x(x), y(y) {}
    bool operator==(const MyPoint& p)const {
		const double eps = 1e-10;
        double d1 = x - p.x;
        double d2 = y - p.y;
        return (((d1 > eps) - (d1 < -eps)) == 0) && (((d2 > eps) - (d2 < -eps)) == 0);
    }
    
;};

class PolyIntersec
{
public:
    int sig(double d) {
        const double eps = 1e-10;
        return(d > eps) - (d < -eps);
    }
    double cross(MyPoint o, MyPoint a, MyPoint b) {
        return(a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
    }
    double area(MyPoint* ps, int n) {
        ps[n] = ps[0];
        double res = 0;
        for (int i = 0; i < n; i++) {
            res += ps[i].x * ps[i + 1].y - ps[i].y * ps[i + 1].x;
        }
        return res / 2.0;
    }
    int lineCross(MyPoint a, MyPoint b, MyPoint c, MyPoint d, MyPoint& p) {
        double s1, s2;
        s1 = cross(a, b, c);
        s2 = cross(a, b, d);
        if (sig(s1) == 0 && sig(s2) == 0) return 2;
        if (sig(s2 - s1) == 0) return 0;
        p.x = (c.x * s2 - d.x * s1) / (s2 - s1);
        p.y = (c.y * s2 - d.y * s1) / (s2 - s1);
        return 1;
    }

    void polygon_cut(MyPoint* p, int& n, MyPoint a, MyPoint b) {
        static MyPoint pp[maxn_poly];
        int m = 0; p[n] = p[0];
        for (int i = 0; i < n; i++) {
            if (sig(cross(a, b, p[i])) > 0) pp[m++] = p[i];
            if (sig(cross(a, b, p[i])) != sig(cross(a, b, p[i + 1])))
                lineCross(a, b, p[i], p[i + 1], pp[m++]);
        }
        n = 0;
        for (int i = 0; i < m; i++)
            if (!i || !(pp[i] == pp[i - 1]))
                p[n++] = pp[i];
        while (n > 1 && p[n - 1] == p[0])n--;
    }

    double intersectArea(MyPoint a, MyPoint b, MyPoint c, MyPoint d) {
        MyPoint o(0, 0);
        int s1 = sig(cross(o, a, b));
        int s2 = sig(cross(o, c, d));
        if (s1 == 0 || s2 == 0)return 0.0;
        if (s1 == -1) swap(a, b);
        if (s2 == -1) swap(c, d);
        MyPoint p[10] = { o,a,b };
        int n = 3;
        polygon_cut(p, n, o, c);
        polygon_cut(p, n, c, d);
        polygon_cut(p, n, d, o);
        double res = fabs(area(p, n));
        if (s1 * s2 == -1) res = -res; return res;
    }

    double intersectArea(MyPoint* ps1, int n1, MyPoint* ps2, int n2) {
        if (area(ps1, n1) < 0) reverse(ps1, ps1 + n1);
        if (area(ps2, n2) < 0) reverse(ps2, ps2 + n2);
        ps1[n1] = ps1[0];
        ps2[n2] = ps2[0];
        double res = 0;
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                res += intersectArea(ps1[i], ps1[i + 1], ps2[j], ps2[j + 1]);
            }
        }
        return res;
    }
};



//MyPoint ps1[maxn], ps2[maxn];
//int n1, n2;
//int main() {
//    PolyIntersec PolyInt;
//    while (scanf("%d%d", &n1, &n2) != EOF) {
//        for (int i = 0; i < n1; i++)
//            scanf("%lf%lf", &ps1[i].x, &ps1[i].y);
//        for (int i = 0; i < n2; i++)
//            scanf("%lf%lf", &ps2[i].x, &ps2[i].y);
//        double ans = PolyInt.intersectArea(ps1, n1, ps2, n2);
//        ans = fabs(PolyInt.area(ps1, n1)) + fabs(PolyInt.area(ps2, n2)) - ans;//ÈÝ³â
//        printf("%.2f\n", ans);
//    }
//    return 0;
//}

