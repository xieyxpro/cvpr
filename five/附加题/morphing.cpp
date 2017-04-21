#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "delaunay.h"
#include "CImg.h"

using namespace std;
using namespace cimg_library;

#define N 3

void split(string& s, string separator, vector<string>& result);
bool inTriangle(PT &a, PT &b, PT &c, PT &p);
void  calAStar(double arcs[N][N], int n, double ans[N][N]);
bool GetMatrixInverse(double src[N][N], int n, double des[N][N]);
double calA(double arcs[N][N], int n);
void multiplyMatrix(double f[N][N], double s[N][N], int n, double ans[N][N]);
int normalizeX(double n);
int normalizeY(double n);
string int2string(int);

int main()
{
	cout << "reading the points\n";
	int x1, y1, x2, y2;
	vector<PT> points1, points2;
	string line;

	ifstream fin("coordinate.txt");
	while (getline(fin, line))
	{
		vector<string> result;
		split(line, " " , result);
		x1 = atoi(result[0].c_str());
		y1 = atoi(result[1].c_str());
		x2 = atoi(result[2].c_str());
		y2 = atoi(result[3].c_str());

		PT p1, p2;
		p1.x = x1;
		p1.y = y1;
		p2.x = x2;
		p2.y = y2;
		points1.push_back(p1);
		points2.push_back(p2);
	}//to get the coordinates of special points
	
	cout << points1.size() << endl;
	cout << "reading finish\n";

	vector<TRIANGLE> t1;
	CMyDelaunay d;
	d.BuildDelaunayEx(points1, t1);//to do delaunay

	cout << t1.size() << endl;

	cout << "finish delaunay\n";

	for (int i = 1; i <= 11; ++i)
	{
		vector<PT> m;

		for (int k = 0; k < t1.size(); ++k)
		{
			PT pm;
			pm.x = i / 12.0 * (points2[t1[k].firstPid-1].x - points1[t1[k].firstPid-1].x) + points1[t1[k].firstPid-1].x;
			pm.y = i / 12.0 * (points2[t1[k].firstPid-1].y - points1[t1[k].firstPid-1].y) + points1[t1[k].firstPid-1].y;
			m.push_back(pm);

			pm.x = i / 12.0 * (points2[t1[k].secendPid-1].x - points1[t1[k].secendPid-1].x) + points1[t1[k].secendPid-1].x;
			pm.y = i / 12.0 * (points2[t1[k].secendPid-1].y - points1[t1[k].secendPid-1].y) + points1[t1[k].secendPid-1].y;
			m.push_back(pm);

			pm.x = i / 12.0 * (points2[t1[k].thirdPid-1].x - points1[t1[k].thirdPid-1].x) + points1[t1[k].thirdPid-1].x;
			pm.y = i / 12.0 * (points2[t1[k].thirdPid-1].y - points1[t1[k].thirdPid-1].y) + points1[t1[k].thirdPid-1].y;
			m.push_back(pm);
		}//cal the coordinate of the middle frame's special points

		cout << m.size() << endl;

		CImg<float> medG(490,700,1,3,0), origin("src1.bmp"), target("src2.bmp");
		for (int u = 0; u < 490; ++u)
		{
			for (int v = 0; v < 700; ++v)
			{
				PT p;
				p.x = u,p.y=v;
				for (int k = 0; k < t1.size(); ++k)
				{
					if (inTriangle(m[k*3], m[k*3+1], m[k*3+2], p)) {
						double m0[N][N] = {
							{m[k*3].x, m[k*3+1].x, m[k*3+2].x},
							{m[k*3].y, m[k*3+1].y, m[k*3+2].y},
							{1, 1, 1}
						};

						double m1[N][N] = {
							{points1[t1[k].firstPid-1].x, points1[t1[k].secendPid-1].x, points1[t1[k].thirdPid-1].x}, 
							{points1[t1[k].firstPid-1].y, points1[t1[k].secendPid-1].y, points1[t1[k].thirdPid-1].y},
							{1, 1, 1}
						};//original picture

						double m2[N][N] = {
							{points2[t1[k].firstPid-1].x, points2[t1[k].secendPid-1].x, points2[t1[k].thirdPid-1].x}, 
							{points2[t1[k].firstPid-1].y, points2[t1[k].secendPid-1].y, points2[t1[k].thirdPid-1].y},
							{1, 1, 1}
						};//target picture

						double m0_inverse[N][N];
						bool flag = GetMatrixInverse(m0, N, m0_inverse);
					    if (false == flag) {
					        cout << "can't get xishu" << endl;
					        return 1;
					    }

					    double tran1[N][N], tran2[N][N];//two transform matrix
					    multiplyMatrix(m1, m0_inverse, N, tran1);
					    multiplyMatrix(m2, m0_inverse, N, tran2);

					    PT p1, p2;
					    p1.x = u * tran1[0][0] + v * tran1[0][1] + tran1[0][2];
					    p1.y = u * tran1[1][0] + v * tran1[1][1] + tran1[1][2];
					    p2.x = u * tran2[0][0] + v * tran2[0][1] + tran2[0][2];
					    p2.y = u * tran2[1][0] + v * tran2[1][1] + tran2[1][2];

					    int px1 = normalizeX(p1.x), py1 = normalizeY(p1.y);//make the coordinate normal
					    int px2 = normalizeX(p2.x), py2 = normalizeY(p2.y);

					    medG(u, v, 0) = (1 - i / 12.0) * origin(px1, py1, 0) + i / 12.0 * target(px2, py2, 0);
					    medG(u, v, 1) = (1 - i / 12.0) * origin(px1, py1, 1) + i / 12.0 * target(px2, py2, 1);
					    medG(u, v, 2) = (1 - i / 12.0) * origin(px1, py1, 2) + i / 12.0 * target(px2, py2, 2);
						break;
					}
				}
			}
		}

		medG.save((int2string(i) + ".bmp").c_str());
	}
	return 0;
}

void split(string& s, string separator, vector<string>& result)
{
	int cutAt;
	while ((cutAt = s.find_first_of(separator)) != s.npos)
	{
		if (cutAt > 0)
		{
			result.push_back(s.substr(0, cutAt));
		}
		s = s.substr(cutAt + 1);
	}
	if (s.length() > 0)
	{
		result.push_back(s);
	}
}

bool inTriangle(PT &a, PT &b, PT &c, PT &p)
{
	double abp = (b.x-a.x)*(p.y-a.y)-(b.y-a.y)*(p.x-a.x);
	double bcp = (c.x-b.x)*(p.y-b.y)-(c.y-b.y)*(p.x-b.x);
	double cap = (a.x-c.x)*(p.y-c.y)-(a.y-c.y)*(p.x-c.x);

	return abp*bcp >= 0 && abp*cap >= 0;
}

bool GetMatrixInverse(double src[N][N], int n, double des[N][N])
{
    double flag = calA(src, n);
    double t[N][N];
    if (0 == flag)
    {
        cout << "can't get the inverse matrix" << endl;
        return false;
    }
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

void multiplyMatrix(double f[N][N], double s[N][N], int n, double ans[N][N])
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			ans[i][j] = 0;
			for (int k = 0; k < n; ++k)
			{
				ans[i][j] += f[i][k] * s[k][j];
			}
		}
	}
}

int normalizeX(double n)
{
	int ans = (int)(n + 0.5);
	if (ans < 0) return 0;
	if (ans >= 490) return 489;
	return ans;
}

int normalizeY(double n)
{
	int ans = (int)(n + 0.5);
	if (ans < 0) return 0;
	if (ans >= 700) return 699;
	return ans;
}

string int2string(int n)
{
	 stringstream stream;  
     stream << n;
     string s;
     stream >> s;
     return s; 
}
