#include "CImg.h"
#include <iostream>
#include <string>
#include <strstream> 
using namespace std;
using namespace cimg_library;

void colorToGray(CImg<unsigned char>& img, CImg<unsigned char>& result);

int main(void)
{
	/*CImg<unsigned char> src;
	src.load_bmp("12.bmp");
	int width = src._width;
	int height = src._height;

	CImg<unsigned char> result(width, height, 1, 1, 0);
	colorToGray(src, result);
	result.save("02.bmp");*/

	
	for (int i = 0; i < 10; ++i)
	{
		CImg<unsigned char> src;
		strstream ss;
		ss << i;
		string s;
		ss >> s;

		string temp = s + "c.bmp";
		string save = s + "g.bmp";
		src.load_bmp(temp.c_str());
		int width = src._width;
		int height = src._height;
		CImg<unsigned char> result(width, height, 1, 1, 0);
		colorToGray(src, result);
		result.save(save.c_str());
	}

	return 0;
}

void colorToGray(CImg<unsigned char>& img, CImg<unsigned char>& result)
{
	cimg_forXY(result, x, y)
	{
		result(x, y) = (img(x, y, 0) * 30 + img(x, y, 1) * 59 + img(x, y, 2) * 11) / 100;
	}
}
