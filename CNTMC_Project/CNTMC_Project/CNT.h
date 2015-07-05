/**
CNT.h
Purpose: Header for CNT.cpp

@author Alex Gabourie
@version 1.00
*/

#include "stdafx.h"
#include <string>

using namespace std;

class CNT
{
	int n; //Hamada n parameter
	int m; //Hamada m parameter
	double length; //Length of entire tube
	double cylinderHeight; //Height of compositional cylinders
	double tubeSeparation; //Separation between compositional cylinders
	double minSpacing; //Minimum spacing from one tube to another
	double diameter; //Diameter of the CNT
	double *positions; //Array of the cylinder and constraint positions

private:
	void setDiameter(int n, int m);

public:
	CNT(string filePath);
	~CNT();
	double getDiameter();
	double getLength();
	double getCylinderHeight();
	double getTubeSeparation();
	double getMinSpacing();
	int getm();
	int getn();

};