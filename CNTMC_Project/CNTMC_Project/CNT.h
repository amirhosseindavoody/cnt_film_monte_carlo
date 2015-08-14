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
#include "segment.h"
#include "chirality.h"


using namespace std;
using namespace Eigen;

struct segment;

class tableElem; //class def to avoid circular dependency

// Extracts and stores all of the pertinent information about a CNT
class CNT
{
	Chirality chir; //Chirality of the nanotube
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
	shared_ptr<vector<shared_ptr<segment>>> calculateSegments(double segLen);
	Vector3d CNT::getPoint(int idx);
	Vector3d CNT::calcEndPt(int idx, double extra);
	Vector3d CNT::calcFinalEndPt(int idx);

public:
	CNT();
	CNT(const string fileName, const string filePath, double segLenMin); //segLenMin in Angstroms
	double getDiameter();
	double getLength();
	double getCylinderHeight();
	double getTubeSeparation();
	double getMinSpacing();
	int getCNTNum();
	bool isInitialized();
	shared_ptr<vector<shared_ptr<segment>>> segs; //The sections of the CNT used to create the MC tables
};

#endif
