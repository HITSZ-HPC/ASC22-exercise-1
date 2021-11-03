#include <iostream>
#include <chrono>
#include "SLIC.h"

typedef std::chrono::high_resolution_clock Clock;

//===========================================================================
/// Load PPM file
///
///
//===========================================================================
void LoadPPM(char *filename, unsigned int **data, int *width, int *height)
{
	char header[1024];
	FILE *fp = NULL;
	int line = 0;

	fp = fopen(filename, "rb");

	// read the image type, such as: P6
	// skip the comment lines
	while (line < 2)
	{
		fgets(header, 1024, fp);
		if (header[0] != '#')
		{
			++line;
		}
	}
	// read width and height
	sscanf(header, "%d %d\n", width, height);

	// read the maximum of pixels
	fgets(header, 20, fp);

	// get rgb data
	unsigned char *rgb = new unsigned char[(*width) * (*height) * 3];
	fread(rgb, (*width) * (*height) * 3, 1, fp);

	*data = new unsigned int[(*width) * (*height) * 4];
	int k = 0;
	for (int i = 0; i < (*height); i++)
	{
		for (int j = 0; j < (*width); j++)
		{
			unsigned char *p = rgb + i * (*width) * 3 + j * 3;
			// a ( skipped )
			(*data)[k] = p[2] << 16; // r
			(*data)[k] |= p[1] << 8; // g
			(*data)[k] |= p[0];		 // b
			k++;
		}
	}

	// ofc, later, you'll have to cleanup
	delete[] rgb;

	fclose(fp);
}

//===========================================================================
/// Load PPM file
///
///
//===========================================================================
int CheckLabelswithPPM(char *filename, int *labels, int width, int height)
{
	char header[1024];
	FILE *fp = NULL;
	int line = 0, ground = 0;

	fp = fopen(filename, "rb");

	// read the image type, such as: P6
	// skip the comment lines
	while (line < 2)
	{
		fgets(header, 1024, fp);
		if (header[0] != '#')
		{
			++line;
		}
	}
	// read width and height
	int w(0);
	int h(0);
	sscanf(header, "%d %d\n", &w, &h);
	if (w != width || h != height)
		return -1;

	// read the maximum of pixels
	fgets(header, 20, fp);

	// get rgb data
	unsigned char *rgb = new unsigned char[(w) * (h)*3];
	fread(rgb, (w) * (h)*3, 1, fp);

	int num = 0, k = 0;
	for (int i = 0; i < (h); i++)
	{
		for (int j = 0; j < (w); j++)
		{
			unsigned char *p = rgb + i * (w)*3 + j * 3;
			// a ( skipped )
			ground = p[2] << 16; // r
			ground |= p[1] << 8; // g
			ground |= p[0];		 // b

			if (ground != labels[k])
				num++;

			k++;
		}
	}

	// ofc, later, you'll have to cleanup
	delete[] rgb;

	fclose(fp);

	return num;
}

//===========================================================================
///	The main function
///
//===========================================================================
int main(int argc, char **argv)
{
	unsigned int *img = NULL;
	int width(0);
	int height(0);

	LoadPPM((char *)"data/case1/input_image.ppm", &img, &width, &height);
	if (width == 0 || height == 0)
		return -1;

	int sz = width * height;
	int *labels = new int[sz];
	int numlabels(0);
	SLIC slic;
	int m_spcount;
	double m_compactness = 10.0;

	// Case 1
	m_spcount = 200;
	auto startTime = Clock::now();
	slic.PerformSLICO_ForGivenK(img, width, height, labels, numlabels, m_spcount, m_compactness); //for a given number K of superpixels
	auto endTime = Clock::now();
	auto compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
	std::cout << "Case 1 Computing time: " << (double)compTime.count() / 1000 << "ms" << std::endl;

	int num = CheckLabelswithPPM((char *)"data/case1/check.ppm", labels, width, height);
	if (num < 0)
	{
		std::cout << "The result for labels is different from output_labels.ppm." << std::endl;
	}
	else
	{
		std::cout << "There are " << num << " points' labels are different from original file." << std::endl;
	}

	// slic.SaveSuperpixelLabels2PPM((char *)"output_labels.ppm", labels, width, height);
	if (labels)
		delete[] labels;
	if (img)
		delete[] img;

	// Case 2
	m_spcount = 400;
	LoadPPM((char *)"data/case2/input_image.ppm", &img, &width, &height);
	if (width == 0 || height == 0)
		return -1;

	sz = width * height;
	labels = new int[sz];

	startTime = Clock::now();
	slic.PerformSLICO_ForGivenK(img, width, height, labels, numlabels, m_spcount, m_compactness); //for a given number K of superpixels
	endTime = Clock::now();
	compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
	std::cout << "Case 2 Computing time: " << (double)compTime.count() / 1000 << "ms" << std::endl;

	num = CheckLabelswithPPM((char *)"data/case2/check.ppm", labels, width, height);
	if (num < 0)
	{
		std::cout << "The result for labels is different from output_labels.ppm." << std::endl;
	}
	else
	{
		std::cout << "There are " << num << " points' labels are different from original file." << std::endl;
	}

	if (labels)
		delete[] labels;
	if (img)
		delete[] img;

	// Case 3
	m_spcount = 150;
	LoadPPM((char *)"data/case3/input_image.ppm", &img, &width, &height);
	if (width == 0 || height == 0)
		return -1;

	sz = width * height;
	labels = new int[sz];

	startTime = Clock::now();
	slic.PerformSLICO_ForGivenK(img, width, height, labels, numlabels, m_spcount, m_compactness); //for a given number K of superpixels
	endTime = Clock::now();
	compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
	std::cout << "Case 3 Computing time: " << (double)compTime.count() / 1000 << "ms" << std::endl;

	num = CheckLabelswithPPM((char *)"data/case3/check.ppm", labels, width, height);
	if (num < 0)
	{
		std::cout << "The result for labels is different from output_labels.ppm." << std::endl;
	}
	else
	{
		std::cout << "There are " << num << " points' labels are different from original file." << std::endl;
	}

	if (labels)
		delete[] labels;
	if (img)
		delete[] img;

	return 0;
}