#ifndef __CNT_H__
#define __CNT_H__

//a = 1.42*sqrt(3) //Amirhossein said ok
#define A_CC 2.459512146747806 //lattice constant CNTs

#include <stdio.h>
#include <string>
#include <vector>
#include <Eigen>
#include <memory>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "tableElem.h"
#include "segment.h"

using namespace std;
using namespace Eigen;

struct segment;

class tableElem; //class def to avoid circular dependency

// Extracts and stores all of the pertinent information about a CNT
class CNT
{
	int n; //Hamada n parameter
	int m; //Hamada m parameter
	double length; //Length of entire tube
	double cylinderHeight; //Height of compositional cylinders
	double tubeSeparation; //Separation between compositional cylinders
	double minSpacing; //Minimum spacing from one tube to another
	double diameter; //Diameter of the CNT
	int cnt_number; //The number associated with the cnt
	vector<vector<double>> positions; //2D array storing positions of cylinders and constraints
	bool initialized = false; //a way to check if variables were initialized
	int number_of_points; // the number of points in the csv file

private:
	void set_diameter(int n, int m);
	void calculate_segments(double segLen);
	void print_segment_points(int first_idx, int n);
	void perform_PCA(int first_idx, int n, double segment_length, vector<double> &first_point, vector<double> &second_point);

public:
	CNT();
	CNT(const string fileName, const string filePath, double segLenMin); //segLenMin in Angstroms
	double getDiameter();
	vector<segment> segments;
};

#endif