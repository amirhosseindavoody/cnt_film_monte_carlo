/**
CNT.h
Purpose: Header for CNT.cpp

@author Alex Gabourie
@version 1.00
*/

#ifndef __CNT_H__
#define __CNT_H__

#include "stdafx.h"
#include <string>
#include <memory>

//a = 1.42*sqrt(3) //Amirhossein said ok
#define A_CC 2.459512146747806 //lattice constant CNTs


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
	int cntNum;
	double **positions; //pointing to pointers
	bool initialized = false; //a way to check if variables were initialized

private:
	void setDiameter(int n, int m);

public:
	CNT();
	CNT(const string fileName,const string filePath);
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
