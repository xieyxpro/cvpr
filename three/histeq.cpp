#include "CImg.h"
#include <iostream>
#include <string>
#include <strstream> 
using std::cout;
using std::endl;
using std::string;
using std::strstream;
using namespace cimg_library;

void grayHisTeq(CImg<unsigned char> &img, CImg<unsigned char> &result);
void colorHisTeq(CImg<unsigned char> &img, CImg<unsigned char> &result);

int main(void)
{
	/*CImg<unsigned char> src;
	src.load_bmp("12.bmp");
	int width = src._width;
	int height = src._height;*/

	/*CImg<unsigned char> result(width, height, 1, 1, 0);
	grayHisTeq(src, result);*/

	/*CImg<unsigned char> result(width, height, 1, 3, 0);
	colorHisTeq(src, result);
	result.save("color.bmp");*/

	/*for (int i = 0; i < 10; ++i)
	{
		CImg<unsigned char> src;
		strstream ss;
		ss << i;
		string s;
		ss >> s;

		string temp = s + "c.bmp";
		string save = s + "c_his.bmp";
		src.load_bmp(temp.c_str());
		int width = src._width;
		int height = src._height;
		CImg<unsigned char> result(width, height, 1, 3, 0);
		colorHisTeq(src, result);
		result.save(save.c_str());
	}*/

	for (int i = 0; i < 10; ++i)
	{
		CImg<unsigned char> src;
		strstream ss;
		ss << i;
		string s;
		ss >> s;

		string temp = s + "g.bmp";
		string save = s + "g_his.bmp";
		src.load_bmp(temp.c_str());
		int width = src._width;
		int height = src._height;
		CImg<unsigned char> result(width, height, 1, 1, 0);
		grayHisTeq(src, result);
		result.save(save.c_str());
	}
	return 0;
}

void grayHisTeq(CImg<unsigned char>& img, CImg<unsigned char>& result)
{
	int width = img._width;
	int height = img._height;
	int total = width * height;

	int gray[256] = {0};
	cimg_forXY(img, x, y)
	{
		if (img(x, y) < 256 && img(x, y) >= 0)
			gray[img(x, y)]++;
	}

	double probability[256];
	for (int i = 0; i < 256; i++)
	{
		probability[i] = gray[i] * 1.0 / total;
	}

	double ccumulate[256];
	ccumulate[0] = probability[0];
	for (int i = 1; i < 256; i++)
	{
		ccumulate[i] = ccumulate[i - 1] + probability[i];
	}

	unsigned char newGrayValue[256];
	for (int i = 0; i < 256; i++)
	{
		newGrayValue[i] = (unsigned char)(255 * ccumulate[i]);
	}

	cimg_forXY(result, x, y)
	{
		result(x, y) = newGrayValue[(int)img(x, y)];
	}
}

void colorHisTeq(CImg<unsigned char>& img, CImg<unsigned char>& result)
{
	CImg<unsigned char> t1(img._width, img._height, 1, 1, 0);
	CImg<unsigned char> t2(img._width, img._height, 1, 1, 0);
	CImg<unsigned char> t3(img._width, img._height, 1, 1, 0);

	cimg_forXY(img, x, y)
	{
		t1(x, y) = img(x, y, 0);
		t2(x, y) = img(x, y, 1);
		t3(x, y) = img(x, y, 2);
	}

	CImg<unsigned char> r1(img._width, img._height, 1, 1, 0);
	CImg<unsigned char> r2(img._width, img._height, 1, 1, 0);
	CImg<unsigned char> r3(img._width, img._height, 1, 1, 0);

	grayHisTeq(t1, r1);
	grayHisTeq(t2, r2);
	grayHisTeq(t3, r3);

	cimg_forXY(result, x, y)
	{
		result(x, y, 0) = r1(x, y);
		result(x, y, 1) = r2(x, y);
		result(x, y, 2) = r3(x, y);
	}
}
