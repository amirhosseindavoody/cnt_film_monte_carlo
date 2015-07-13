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
#include "exciton.h"
#include <iostream>

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
	shared_ptr<exciton> ex1; //first energy level for exciton
	shared_ptr<exciton> ex2; //second energy level for exciton

	/**
	Determines if the segment has an exciton of the same type
	as the passed exciton

	@param e The exciton desired to see if a slot is available
	@return Whether or not the exciton can be added
	*/
	bool hasExciton(shared_ptr<exciton> e)
	{
		int currEn = e->getEnergy();
		//parameter check
		if (currEn != 1 || currEn != 2)
		{
			cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}

		//if energy is 1 and that slot is occupied return true
		if (currEn == 1 && !(ex1 == nullptr))
		{
			return true;
		}
		//if energy is 2 and that slot is occupied return true
		if (currEn == 2 && !(ex2 == nullptr))
		{
			return true;
		}
		return false; //slot is empty
	}

	/**
	Sets the exciton in the correct slot to the exciton passed to function.
	

	@param e The exciton desired to be added to the segment
	@return True if assignment works and false if assignment not successful
	*/
	bool setExciton(shared_ptr<exciton> e)
	{
		int currEn = e->getEnergy();
		//parameter check
		if (currEn != 1 || currEn != 2)
		{
			cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}

		if (this->hasExciton(e))
		{
			if (currEn == 1)
			{
				ex1 = e;
				return true;
			}
			if (currEn == 2)
			{
				ex2 = e;
				return true;
			}
		}
		return false;
	}

	/**
	Removes the exciton of the correct type 
	
	@param e The exciton desired to be removed from the segment
	@return True if exciton removed, false if no exciton to remove
	*/
	bool removeExciton(shared_ptr<exciton> e)
	{
		int currEn = e->getEnergy();
		//parameter check
		if (currEn != 1 || currEn != 2)
		{
			cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}

		if (this->hasExciton(e))
		{
			if (currEn == 1)
			{
				ex1 = nullptr;
				return true;
			}
			if (currEn == 2)
			{
				ex2 = nullptr;
				return true;
			}
		}
		return false; //no exciton present to remove
	}

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
	Vector3d CNT::calcFinalEndPt(int idx);

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
