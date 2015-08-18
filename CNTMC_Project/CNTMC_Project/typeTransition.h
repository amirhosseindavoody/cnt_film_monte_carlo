/**
typeTransition.h
Purpose: header for typeTransition.cpp

@author Alex Gabourie
@version 1.00
*/

#pragma once
#include "chirality.h"
#include <memory>
#include <vector>
#include "energyTransition.h"

using namespace std;

class typeTransition
{

	Chirality source_chirality; //chirality of tube exciton is currently located
	Chirality dest_chirality; //chirality of tube exciton will move to
	shared_ptr<vector<energyTransition>> ETransList; //List of the posible energy transitions

public:
	typeTransition();
	typeTransition(Chirality &src, Chirality &dest, uint32_t r_length, uint32_t theta_length);
	~typeTransition();
	shared_ptr<vector<energyTransition>> getETransList();
};

