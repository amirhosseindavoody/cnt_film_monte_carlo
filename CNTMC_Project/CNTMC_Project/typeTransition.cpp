/**
typeTransition.cpp
Purpose: Store transition table data for a particular chirality to another 
particular chirality.

@author Alex Gabourie
@version 1.00
*/

#include "stdafx.h"
#include "typeTransition.h"
#include "energyTransition.h"
#include <array>

using namespace std;

typeTransition::typeTransition()
{
	ETransList = make_shared<vector<energyTransition>>(vector<energyTransition>(0));
}

typeTransition::typeTransition(Chirality &src, Chirality &dest, uint32_t r_length, uint32_t theta_length)
{
	source_chirality = src;
	dest_chirality = dest;
	ETransList = make_shared<vector<energyTransition>>(vector<energyTransition>(0));
	//initialize array
	std::array<

}


typeTransition::~typeTransition()
{
}
