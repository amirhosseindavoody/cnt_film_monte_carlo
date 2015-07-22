/**
exciton.cpp
Purpose: Stores relevant exciton information

@author Alex Gabourie
@version 1.00
*/

#include "stdafx.h"
#include "exciton.h"

/**
Creates exciton object

@return tableElem Object
*/
exciton::exciton()
{
	
}

/**
Creates exciton object

@param cidx Index of the CNT the exciton belongs to
@param sidx Index of the segment the exciton belongs to
@param energy Whether the 1st or 2nd energy level
@return tableElem Object
*/
exciton::exciton(int cidx, int sidx, int energy)
{
	cntidx = cidx;
	segidx = sidx;
	energyNum = energy;
}

/**
destroys exciton object

@return tableElem Object
*/
exciton::~exciton()
{
}

/**
Sets cnt index

@param cidx Index of the CNT the exciton belongs to
*/
void exciton::setCNTidx(int cidx)
{
	cntidx = cidx;
}

/**
Sets segment index

@param sidx Index of the segment the exciton belongs to
*/
void exciton::setSegidx(int sidx)
{
	segidx = sidx;
}

/**
Sets energy level

@param energy Whether the 1st or 2nd energy level
*/
void exciton::setEnergy(int energy)
{
	energyNum = energy;
}

/**
Sets boolean for whether or not at end contact

@param atContact Whether or not at end contact
*/
void exciton::setAtOutContact(bool atContact)
{
	atOutContact = atContact;
}


/**
Gets cnt index

@return cnt index
*/
int exciton::getCNTidx()
{
	return cntidx;
}

/**
Get segment index

@return segment index
*/
int exciton::getSegidx()
{
	return segidx;
}

/**
Gets energy level

@return energy level
*/
int exciton::getEnergy()
{
	return energyNum;
}

/**
Checks whether or not the exciton is at the exit contact
*/
bool exciton::isAtOutContact()
{
	return atOutContact;
}
