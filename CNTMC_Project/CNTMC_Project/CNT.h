/**
CNT.h
Purpose: Header for CNT.cpp

@author Alex Gabourie
@version 1.00
*/

#ifndef __CNT_H__
#define __CNT_H__

//a = 1.42*sqrt(3) //Amirhossein said ok
#define A_CC 2.459512146747806 //lattice constant CNTs

#include "stdafx.h"
#include <string>
#include <vector>
#include <Eigen>
#include <memory>
#include "tableElem.h"

using namespace std;
using namespace Eigen;

class tableElem; //class def to avoid circular dependency

//Stores the position information about 
struct segment
{
	int segNum;
	Vector3d p1; //first point in segment
	Vector3d p2; //second point in segment
	Vector3d mid; //middle point in segment
	shared_ptr<vector<tableElem>> tbl;
};

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
	int cntNum; //The number associated with the cnt
	vector<vector<double>> positions; //2D array storing positions of cylinders and constraints
	bool initialized = false; //a way to check if variables were initialized
	int numPt; // the number of points in the csv file

private:
	void setDiameter(int n, int m);
	shared_ptr<vector<segment>> calculateSegments(double segLen);
	Vector3d CNT::getPoint(int idx);
	Vector3d CNT::calcEndPt(int idx, double extra);

public:
	CNT();
	CNT(const string fileName, const string filePath, double segLenMin); //segLenMin in Angstroms
	double getDiameter();
	double getLength();
	double getCylinderHeight();
	double getTubeSeparation();
	double getMinSpacing();
	int getm();
	int getn();
	int getCNTNum();
	bool isInitialized();
	shared_ptr<vector<segment>> segs; //The sections of the CNT used to create the MC tables
};

#endif
