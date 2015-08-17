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