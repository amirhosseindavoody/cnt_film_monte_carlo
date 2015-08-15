/**
exciton.h
Purpose: Header for exciton.cpp

@author Alex Gabourie
@version 1.00
*/
#pragma once
#include "energyTransition.h"

class exciton
{
	int cntidx; //Index of the CNT the exciton belongs to
	int segidx; //Index of the segment the exciton belongs to
	bool atOutContact; //States whether or not the exciton is ready to leave mesh
	double textra; //The extra amount of time past deltaT the exciton had in last step
	energy e;

public:
	exciton();
	exciton(int cidx, int sidx, energy energy_new);
	~exciton();
	void setCNTidx(int cidx);
	void setSegidx(int sidx);
	void setEnergy(energy e_new);
	void setAtOutContact(bool atContact);
	void setTExtra(double t);
	int getCNTidx();
	int getSegidx();
	energy getEnergy();
	bool isAtOutContact();
	double getTExtra();

};

