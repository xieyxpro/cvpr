#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;
const double EP = 0.00000001;

// 点结构
typedef struct PT{
    double x;
    double y;
    PT(){}
    PT(double a, double b){x=a,y=b;}
}PT;


typedef struct TRIANGLE
{
    int firstPid, secendPid, thirdPid;
    TRIANGLE();
    TRIANGLE(int f, int s, int t){firstPid=f;secendPid=s;thirdPid=t;}
}TRIANGLE;

class CMyDelaunay  
{
public:
    static void BuildDelaunayEx(vector<PT>& vecPT, vector<TRIANGLE>& vecTriangleWork);
public:
    CMyDelaunay() {}
    ~CMyDelaunay() {}
};

void CMyDelaunay::BuildDelaunayEx(vector<PT>& vecPT, vector<TRIANGLE>& vecTriangleWork)
{
    cout << "start delaunaying\n";
    ifstream fin("delaunay_result.txt", ios::in);
    int f, s, t;
    while (fin >> f >> s >> t)
    {
        vecTriangleWork.push_back(TRIANGLE(f,s,t));
    }
    cout << "finish delaunaying\n";
}
