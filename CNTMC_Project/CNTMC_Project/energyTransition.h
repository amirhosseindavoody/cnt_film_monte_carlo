/**
energyTransition.h
Purpose: header for energyTransition.cpp

@author Alex Gabourie
@version 1.00
*/

#pragma once
#include <memory>
#include <vector>
#include "point.h"

using namespace std;

//Type to store energy information
enum energy
{
	E11,
	E22
};

class energyTransition
{

	energy source_energy; //energy of exciton on source segment
	energy dest_energy; //energy of exciton on destination segment
	shared_ptr<vector<vector<double>>> rtTable; // r and theta table values
	uint32_t r_len; //number of r values
	uint32_t t_len; //number of theta values

public:
	energyTransition();
	energyTransition(energy src, energy dest, uint32_t r_length, uint32_t theta_length);
	~energyTransition();
	void initializeTable(uint32_t r_length, uint32_t theta_length);
	void setDestEnergy(energy dest);
	void setSrcEnergy(energy src);
	energy getDestEnergy();
	energy getSrcEnergy();
	void setTableValue(uint32_t r_idx, uint32_t t_idx, double val);
	double getTableValue(uint32_t r_idx, uint32_t t_idx);
	void interpolateRate(point &p00, point &p01, point &p10, point &p11, point &interp);
};

