﻿/**
segment.h
Purpose: Header for segment.cpp

@author Alex Gabourie
@version 1.00
*/
#pragma once

#include "tableElem.h"
#include "exciton.h"
#include <Eigen>
#include <memory>

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
	shared_ptr<exciton> ex1; //first energy level for exciton
	shared_ptr<exciton> ex2; //second energy level for exciton

	//instance methods
	bool hasExciton(shared_ptr<exciton> e);
	bool setExciton(shared_ptr<exciton> e);
	bool removeExciton(shared_ptr<exciton> e);
};