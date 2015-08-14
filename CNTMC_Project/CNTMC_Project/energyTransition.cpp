/**
energyTransition.cpp
Purpose: Stores energies that exciton will transfer from and to 
as well as the rate based on the r, theta tables

@author Alex Gabourie
@version 1.00
*/

#include "stdafx.h"
#include "energyTransition.h"

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
