#include <iostream>
#include <fstream>
#include "CImg.h"
#include <cmath>
using namespace std;
using namespace cimg_library;

int main(void)
{
	CImg<float> src1("src1.bmp"), src2("src2.bmp");
	ofstream fout1("src1.txt"), fout2("src2.txt");

	cimg_forXY(src1, x, y)
	{
		fout1 << src1(x,y,0) << " " << src1(x,y,1) << " " << src1(x,y,2) << endl;
		fout2 << src2(x,y,0) << " " << src2(x,y,1) << " " << src2(x,y,2) << endl;
	}
	return 0;
}

