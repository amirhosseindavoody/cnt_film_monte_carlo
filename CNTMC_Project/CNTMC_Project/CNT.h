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
#include <array>
#include <vector>

using namespace std;

//Stores the position information about 
struct segment
{
	vector<double> p1 = vector<double>(3); //beginning of section point
	vector<double> p2 = vector<double>(3); //end of section point
	vector<double> mid = vector<double>(3); //middle of section point
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
	array<vector<double>,3> positions; //2D array storing positions of cylinders and constraints
	bool initialized = false; //a way to check if variables were initialized
	vector<segment> segs; //The sections of the CNT used to create the MC tables

private:
	void setDiameter(int n, int m);
	vector<segment> calculateSegments(double segLen);

public:
	CNT();
	CNT(const string fileName,const string filePath, double segLen);
	double getDiameter();
	double getLength();
	double getCylinderHeight();
	double getTubeSeparation();
	double getMinSpacing();
	int getm();
	int getn();
	int getCNTNum();
	bool isInitialized();

};



#endif
