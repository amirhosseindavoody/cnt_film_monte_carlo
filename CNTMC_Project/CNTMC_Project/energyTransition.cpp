/**
energyTransition.cpp
Purpose: Stores energies that exciton will transfer from and to 
as well as the rate based on the r, theta tables

@author Alex Gabourie
@version 1.00
*/

#include "stdafx.h"
#include "energyTransition.h"
#include "point.h"

/**
Sets the energyTransition object to some default values. DO NOT USE CONSTRUCTED
THIS WAY. This is only to appease the compiling gods.

@return energyTransition Object
*/
energyTransition::energyTransition()
{
	//default values chosen. table must be initialized
}

/**
Constructs energy transition object and sets up the r theta table 
*/
energyTransition::energyTransition(energy src,energy dest, uint32_t r_length, uint32_t theta_length)
{
	//init energies
	source_energy = src;
	dest_energy = dest;
	//init table
	r_len = r_length;
	t_len = theta_length;
	rtTable = make_shared<vector<vector<double>>>(vector<vector<double>>(r_length));
	//set table to correct size
	for (auto it = rtTable->begin(); it != rtTable->end(); ++it)
	{
		it->resize(theta_length);
	}
}



/**
Destructor for energyTransition object
*/
energyTransition::~energyTransition()
{
}

/**
Initializes the r theta table. Removes any old data stored in table

@param r_length The number of r values
@param theta_length The number of theta values
*/
void energyTransition::initializeTable(uint32_t r_length, uint32_t theta_length)
{
	//init table
	rtTable = make_shared<vector<vector<double>>>(vector<vector<double>>(r_length));
	//set table to correct size
	r_len = r_length;
	t_len = theta_length;
	for (auto it = rtTable->begin(); it != rtTable->end(); ++it)
	{
		it->resize(theta_length);
	}
}

/**
Returns the destination energy

@return Destination energy
*/
energy energyTransition::getDestEnergy()
{
	return dest_energy;
}

/**
Returns the source energy

@return source energy
*/
energy energyTransition::getSrcEnergy()
{
	return source_energy;
}

/**
Sets the destination energy

@param Destination energy
*/
void energyTransition::setDestEnergy(energy dest)
{
	dest_energy = dest;
}

/**
Sets the source energy

@param source energy
*/
void energyTransition::setSrcEnergy(energy src)
{
	source_energy = src;
}

/**
Sets the given index of the table to the value val

@param r_idx The index of the r value
@param t_idx The index of the theta value
@param val The value to place at the r_idx and t_idx to change the rate
*/
void energyTransition::setTableValue(uint32_t r_idx, uint32_t t_idx, double val)
{
	if (r_idx >= r_len && t_idx >= t_len)
	{
		printf("Invalid index for r and theta table.\n");
		system("pause");
		exit(EXIT_FAILURE);
	} 
	else if (val < 0)
	{
		printf("Invalid transition rate passed to r and theta table. " 
			"Transition rate must be positive.\n");
		system("pause");
		exit(EXIT_FAILURE);
	}
	//set value
	(*rtTable)[r_idx][t_idx] = val;
}

/**
Gets the transition rate value from the r theta table at the specified index

@param r_idx The index of the r value
@param t_idx The index of the theta value
@return value at the input index
*/
double energyTransition::getTableValue(uint32_t r_idx, uint32_t t_idx)
{
	if (r_idx >= r_len && t_idx >= t_len)
	{
		printf("Invalid index for r and theta table.\n");
		system("pause");
		exit(EXIT_FAILURE);
	}
	//return rate value
	return (*rtTable)[r_idx][t_idx];
}

/**
R and theta values coming from two segments will most likely not match an r/theta
pair from the rate tables. In this case a rate must be interpolated.


p00 ----------- p01
 |               |
 |               |
 |               |
 |               |
p10 ----------- p11 

Current algorithm is bilinear interpolation

@param p00
@param p01
@param p10
@param p11
@param interp
*/
void energyTransition::interpolateRate(point &p00, point &p01, point &p10, point &p11, point &interp)
{
	//get point values
	p00.val = getTableValue(p00.r_idx, p00.t_idx);
	p01.val = getTableValue(p01.r_idx, p01.t_idx);
	p10.val = getTableValue(p10.r_idx, p10.t_idx);
	p11.val = getTableValue(p11.r_idx, p11.t_idx);

	auto r1 = p00.r_val;
	auto r2 = p11.r_val;
	auto t1 = p00.t_val;
	auto t2 = p11.t_val;
	auto r = interp.r_val;
	auto t = interp.t_val;

	interp.val = 1 / ((r2 - r1)*(t2 - t1))*(p00.val*(r2 - r)*(t2 - t) +
		p10.val*(r - r1)*(t2 - t) + p01.val*(r2 - r)*(t - t1) + p11.val*(r - r1)*(t - t1));
}

