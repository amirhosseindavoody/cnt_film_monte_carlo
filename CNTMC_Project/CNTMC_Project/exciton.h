/**
exciton.h
Purpose: Header for exciton.cpp

@author Alex Gabourie
@version 1.00
*/
#pragma once
class exciton
{
	int cntidx; //Index of the CNT the exciton belongs to
	int segidx; //Index of the segment the exciton belongs to
	int energyNum; //Whether the 1st or 2nd energy level
	bool atOutContact; //States whether or not the exciton is ready to leave mesh
	double textra; //The extra amount of time past deltaT the exciton had in last step
	bool justInjected;

public:
	exciton();
	exciton(int cidx, int sidx, int energy);
	~exciton();
	void setCNTidx(int cidx);
	void setSegidx(int sidx);
	void setEnergy(int energy);
	void setAtOutContact(bool atContact);
	void setTExtra(double t);
	int getCNTidx();
	int getSegidx();
	int getEnergy();
	bool isAtOutContact();
	double getTExtra();

};

