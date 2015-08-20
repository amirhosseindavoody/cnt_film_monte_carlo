/**
paramStructs.h
Purpose: Stores relevant information for updating the rate tables for
both the calculating and table reading methods

@author Alex Gabourie
@version 1.00
*/

#pragma once
#include <vector>
#include <memory>
#include "typeTransition.h"
#include "segment.h"

using namespace std;

/*
seg The segment the table element is being added to
r The separation between two segs
theta Angle between two segs
i CNT index
j segment index
rate_tot The running sum of the rates for the segment
c2c Amirhossein's tables used to extract rates
r_vec The vector of r's used to calculate Amirhossein's table index
t_vec The vector of thetas used to caluclate Amirhossein's table index
CNT_List The list of CNTs that stores all segments and rate tables
*/

struct tableUpdater
{
	vector<shared_ptr<segment>>::iterator seg;
	double r;
	double theta;
	int i;
	int j;
	double rate_tot;
	shared_ptr<vector<vector<typeTransition>>> c2c;
	shared_ptr<vector<double>> r_vec;
	shared_ptr<vector<double>> t_vec;
	shared_ptr<vector<CNT>> CNT_List;
};

struct heatMapInfo
{
	shared_ptr<vector<vector<int>>> map;
	shared_ptr<vector<double>> rs;
	shared_ptr<vector<double>> thetas;
};

struct chirPair
{
	Chirality c1;
	Chirality c2;
};

