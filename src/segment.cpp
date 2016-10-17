/**
segment.cpp
Purpose: Segment struct used in each CNT object

@author Alex Gabourie
@version 1.00
*/
#include <stdio.h>
#include <iostream>

#include "segment.h"

/**
Determines if the segment has an exciton of the same type
as the passed exciton

@param e The exciton desired to see if a slot is available
@return Whether or not the exciton can be added
*/
bool segment::hasExciton(shared_ptr<exciton> e)
{
	int currEn = e->getEnergy();
	//parameter check
	if (currEn != 1 && currEn != 2)
	{
		cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	//if energy is 1 and that slot is occupied return true
	if (currEn == 1 && !(ex1 == nullptr))
	{
		return true;
	}
	//if energy is 2 and that slot is occupied return true
	if (currEn == 2 && !(ex2 == nullptr))
	{
		return true;
	}
	return false; //slot is empty
}

/**
Sets the exciton in the correct slot to the exciton passed to function.


@param e The exciton desired to be added to the segment
@return True if assignment works and false if assignment not successful
*/
bool segment::setExciton(shared_ptr<exciton> e)
{
	int currEn = e->getEnergy();
	//parameter check
	if (currEn != 1 && currEn != 2)
	{
		cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	if (!this->hasExciton(e))
	{
		if (currEn == 1)
		{
			ex1 = e;
			return true;
		}
		if (currEn == 2)
		{
			ex2 = e;
			return true;
		}
	}
	return false;
}

/**
Removes the exciton of the correct type

@param e The exciton desired to be removed from the segment
@return True if exciton removed, false if no exciton to remove
*/
bool segment::removeExciton(shared_ptr<exciton> e)
{
	int currEn = e->getEnergy();
	//parameter check
	if (currEn != 1 && currEn != 2)
	{
		cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	if (this->hasExciton(e))
	{
		if (currEn == 1)
		{
			ex1 = nullptr;
			return true;
		}
		if (currEn == 2)
		{
			ex2 = nullptr;
			return true;
		}
	}
	return false; //no exciton present to remove
}

/**
Checks to see if the exciton that is passes is the exact exciton that
already exists in the location.
*/
bool segment::hasExactExciton(shared_ptr<exciton> e)
{
	//If adding the same exciton to the same location, that is self scattering and allowed
	if (e == ex1 || e == ex2)
	{
		return true;
	}
	return false;
}
