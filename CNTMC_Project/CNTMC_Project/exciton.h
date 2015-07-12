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

public:
	exciton();
	exciton(int cidx, int sidx, int energy);
	~exciton();
	void setCNTidx(int cidx);
	void setSegidx(int sidx);
	void setEnergy(int energy);
	int getCNTidx();
	int getSegidx();
	int getEnergy();

};

