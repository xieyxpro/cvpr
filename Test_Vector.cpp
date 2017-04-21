// Test_Vector.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"

#include <vector>
#include <algorithm>

#include <iostream>
#include "CImg.h"
using namespace cimg_library;
using namespace std;

void test_cimg()
{
	CImg<unsigned char> SrcImg;
	SrcImg.load_bmp("1.bmp");

	SrcImg.display();//first: to display the image
	
	int w = SrcImg._width;
	int h = SrcImg._height;
	
	CImg<unsigned char> TempImg1(w, h, 1, 3, 0), TempImg2(w, h, 1, 3, 0), TempImg3(w, h, 1, 3, 0);
	
	cimg_forXY(SrcImg, x, y)
	{
		if (SrcImg(x, y, 0) == 255 && SrcImg(x, y, 1) == 255 && SrcImg(x, y, 2) == 255)
		{	//turn white to red
			TempImg1(x, y, 0) = 255;
			TempImg1(x, y, 1) = 0;
			TempImg1(x, y, 2) = 0;
		}
		else if (SrcImg(x, y, 0) == 0 && SrcImg(x, y, 1) == 0 && SrcImg(x, y, 2) == 0)
		{	//turn black to green
			TempImg1(x, y, 0) = 0;
			TempImg1(x, y, 1) = 255;
			TempImg1(x, y, 2) = 0;
		}
		else {
			TempImg1(x, y, 0) = SrcImg(x, y, 0);
			TempImg1(x, y, 1) = SrcImg(x, y, 1);
			TempImg1(x, y, 2) = SrcImg(x, y, 2);
		}
	}

	TempImg1.save("2.bmp");//second

	for (int i = 20; i <= 80; ++i)
	{
		for (int j = 20; j <= 80; ++j)
		{
			double dis = sqrt((double)((i - 50) * (i - 50) + (j - 50) * (j - 50)));
			if (dis <= 30.0)
			{
				TempImg2(i, j, 0) = 0;
				TempImg2(i, j, 1) = 0;
				TempImg2(i, j, 2) = 255;
			}
		}
	}
	TempImg2.save("3.bmp");//third

	for (int i = 47; i <= 53; ++i)
	{
		for (int j = 47; j <= 53; ++j)
		{
			double dis = sqrt((double)((i - 50) * (i - 50) + (j - 50) * (j - 50)));
			if (dis <= 3.0)
			{
				TempImg3(i, j, 0) = 255;
				TempImg3(i, j, 1) = 255;
				TempImg3(i, j, 2) = 0;
			}
		}
	}
	TempImg3.save("4.bmp");//fourth
}

int main(void)
{
	test_cimg();

	return 0;
}

