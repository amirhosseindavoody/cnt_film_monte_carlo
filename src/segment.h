/**
segment.h
Purpose: Header for segment.cpp

@author Alex Gabourie
@version 1.00
*/
#pragma once

#include <Eigen>
#include <memory>

#include "tableElem.h"
#include "exciton.h"

using namespace std;
using namespace Eigen;

class tableElem;

//Stores the position information about 
struct segment
{
	//instance variables
	int segNum;
	Vector3d p1; //first point in segment
	Vector3d p2; //second point in segment
	Vector3d mid; //middle point in segment
	shared_ptr<vector<tableElem>> tbl;
	shared_ptr<vector<double>> rateVec; // Vector to store all rates
	
	/*
	7/20/15: It was decided that there are no limits on the number of excitons that
	can be on a segment. Exciton structures and instance methods are left in for the
	case that excitons are limited to one of each type for each segment
	*/
	shared_ptr<exciton> ex1; //first energy level for exciton
	shared_ptr<exciton> ex2; //second energy level for exciton

	//instance methods
	bool hasExciton(shared_ptr<exciton> e);
	bool setExciton(shared_ptr<exciton> e);
	bool removeExciton(shared_ptr<exciton> e);
	bool hasExactExciton(shared_ptr<exciton> e);
};
