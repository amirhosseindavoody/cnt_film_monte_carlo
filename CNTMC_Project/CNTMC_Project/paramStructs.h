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
#include <array>

using namespace std;


struct tableUpdater
{
	vector<shared_ptr<segment>>::iterator seg; //The segment the table is being added to
	double r; // The separation between the two segs
	double theta; //Angle between two segs
	int src_cnt; //source CNT Index
	int src_seg; //source segment index
	int dest_cnt;// destination CNT index
	int dest_seg; //destionation segment index
	array<double, 2> rate_tot; //running sum of the rates for both source energies
	shared_ptr<vector<vector<typeTransition>>> c2c; //AH tables used to extract rates
	shared_ptr<vector<double>> r_vec; //The vector of r's used to calculate Amirhossein's table index
	shared_ptr<vector<double>> t_vec; //The vector of thetas used to caluclate Amirhossein's table index
	shared_ptr<vector<CNT>> CNT_List; //The list of CNTs that stores all segments and rate tables
	shared_ptr<vector<Chirality>> chirList; //List of chiralities in simulation
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

