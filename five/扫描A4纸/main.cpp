#include <iostream>
#include <vector>
#include <cmath>
#include "CImg.h"
using namespace std;
using namespace cimg_library;

#define N 8

void  calAStar(double arcs[N][N], int n, double ans[N][N]);
bool GetMatrixInverse(double src[N][N], int n, double des[N][N]);
double calA(double arcs[N][N], int n);//calculate det

int main(void)
{
	vector<pair<int, int> > P, newPoint;
	P.push_back(make_pair(0, 0));
	P.push_back(make_pair(299, 0));
	P.push_back(make_pair(0, 419));
	P.push_back(make_pair(299, 419));

    cout << "please input four points\n";
	for (int i = 0; i < 4; ++i)
	{
		int x, y;
		cin >> x >> y;
		newPoint.push_back(make_pair(x,y));
	}

	int uv[8] = { newPoint[0].first, newPoint[0].second,
        newPoint[1].first, newPoint[1].second,
        newPoint[2].first, newPoint[2].second,
        newPoint[3].first, newPoint[3].second };
    
    double src[8][8] =
    { { P[0].first, P[0].second, 1, 0, 0, 0, -newPoint[0].first*P[0].first, -newPoint[0].first*P[0].second },
    { 0, 0, 0, P[0].first, P[0].second, 1, -newPoint[0].second*P[0].first, -newPoint[0].second*P[0].second },

    { P[1].first, P[1].second, 1, 0, 0, 0, -newPoint[1].first*P[1].first, -newPoint[1].first*P[1].second },
    { 0, 0, 0, P[1].first, P[1].second, 1, -newPoint[1].second*P[1].first, -newPoint[1].second*P[1].second },

    { P[2].first, P[2].second, 1, 0, 0, 0, -newPoint[2].first*P[2].first, -newPoint[2].first*P[2].second },
    { 0, 0, 0, P[2].first, P[2].second, 1, -newPoint[2].second*P[2].first, -newPoint[2].second*P[2].second },

    { P[3].first, P[3].second, 1, 0, 0, 0, -newPoint[3].first*P[3].first, -newPoint[3].first*P[3].second },
    { 0, 0, 0, P[3].first, P[3].second, 1, -newPoint[3].second*P[3].first, -newPoint[3].second*P[3].second } };

	double matrix_after[N][N]{};
    bool flag = GetMatrixInverse(src, N, matrix_after);
    if (false == flag) {
        cout << "can't get xishu" << endl;
        return 1;
    }

    double xs[8];
    for (int i = 0; i < 8; i++) {
        double sum = 0;
        for (int t = 0; t < 8; t++) {
            sum += matrix_after[i][t] * uv[t];
        }
        xs[i] = sum;
    }

    for (int i = 0; i<3; i++)
    {
        for (int j = 0; j<3; j++)
        {
            cout << xs[i * 3 + j] << " ";
        }
        cout << endl;
    }

    CImg<float> paint("test.bmp"), outputimg(300,420,1,3,0);
    cimg_forXY(outputimg, x, y) {
        double px = xs[0] * x + xs[1] * y + xs[2];
        double py = xs[3] * x + xs[4] * y + xs[5];
        double p = xs[6] * x + xs[7] * y + 1;

        double u = px / p;
        double v = py / p;

        int uu, vv;
        if (floor(u) < 0) uu = 0;
        else uu = floor(u);

        if (u + 1 > 299) uu = 298;

        if (floor(v) < 0) vv = 0;
        else vv = floor(v);

        if (v + 1 > 419) vv = 418;

        double t = u-uu, s = v-vv;
        outputimg(x,y,0)=(1-t)*(1-s)*paint(uu,vv,0)+(1-t)*s*paint(uu,vv+1,0)+t*(1-s)*paint(uu+1,vv,0)+t*s*paint(uu+1,vv+1,0);
        outputimg(x,y,1)=(1-t)*(1-s)*paint(uu,vv,1)+(1-t)*s*paint(uu,vv+1,1)+t*(1-s)*paint(uu+1,vv,1)+t*s*paint(uu+1,vv+1,1);
        outputimg(x,y,2)=(1-t)*(1-s)*paint(uu,vv,2)+(1-t)*s*paint(uu,vv+1,2)+t*(1-s)*paint(uu+1,vv,2)+t*s*paint(uu+1,vv+1,2);
        /*outputimg(x, y, 0) = paint(u, v, 0);
        outputimg(x, y, 1) = paint(u, v, 1);
        outputimg(x, y, 2) = paint(u, v, 2);*/
    }

    outputimg.save("output.bmp");
    system("pause");
    
	return 0;
}

bool GetMatrixInverse(double src[N][N], int n, double des[N][N])
{
    double flag = calA(src, n);
    double t[N][N];
    if (0 == flag)
    {
        cout << "can't get the inverse matrix" << endl;
        return false;
    else
    {
        calAStar(src, n, t);
        for (int i = 0; i<n; i++)
        {
            for (int j = 0; j<n; j++)
            {
                des[i][j] = t[i][j] / flag;
            }

        }
    }

    return true;
}

void  calAStar(double arcs[N][N], int n, double ans[N][N])
{
    if (n == 1)
    {
        ans[0][0] = 1;
        return;
    }
    int i, j, k, t;
    double temp[N][N];
    for (i = 0; i<n; i++)
    {
        for (j = 0; j<n; j++)
        {
            for (k = 0; k<n - 1; k++)
            {
                for (t = 0; t<n - 1; t++)
                {
                    temp[k][t] = arcs[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
                }
            }


            ans[j][i] = calA(temp, n - 1);
            if ((i + j) % 2 == 1)
            {
                ans[j][i] = -ans[j][i];
            }
        }
    }
}

double calA(double arcs[N][N], int n)
{
    if (n == 1)
    {
        return arcs[0][0];
    }
    double ans = 0;
    double temp[N][N] = { 0.0 };
    int i, j, k;
    for (i = 0; i<n; i++)
    {
        for (j = 0; j<n - 1; j++)
        {
            for (k = 0; k<n - 1; k++)
            {
                temp[j][k] = arcs[j + 1][(k >= i) ? k + 1 : k];

            }
        }
        double t = calA(temp, n - 1);
        if (i % 2 == 0)
        {
            ans += arcs[0][i] * t;
        }
        else
        {
            ans -= arcs[0][i] * t;
        }
    }
    return ans;
}
