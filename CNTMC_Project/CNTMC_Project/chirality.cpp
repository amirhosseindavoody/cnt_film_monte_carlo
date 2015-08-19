/**
chirality.cpp
Purpose: Stores the hamada parameters n, m as an object.
This allows for better organization when dealing with 
preparing tables based on exciton transition rate data.

@author Alex Gabourie
@version 1.00
*/
#include "stdafx.h"
#include "chirality.h"
#include <exception>


/**
Sets the Chirality object to some default values. DO NOT USE Chirality CONSTRUCTED
THIS WAY. This is only to appease the compiling gods.

@return Chirality Object
*/
Chirality::Chirality()
{
	n = 0;
	m = 0;
}

/**
Destructor for chirality object
*/
Chirality::~Chirality() {}

/**
Gets m parameter of the CNT

@param void
@return m
*/
int Chirality::getm() const
{
	return m;
}

/**
Gets n parameter of the CNT

@param void
@return n
*/
int Chirality::getn() const
{
	return n;
}


/**
Sets the n parameter of the chirality object

@param n The hamada n parameter
*/
void Chirality::setn(int n_new)
{
	n = n_new;
}


/**
Sets the m parameter of the chirality object

@param m The hamada n parameter
*/
void Chirality::setm(int m_new)
{
	if (m_new > n)
	{
		printf("Chirality Error: m must be less than or equal to n.\n");
		system("pause");
		exit(EXIT_FAILURE);
	}
	m = m_new;
}

/**
Compares this object to the passed object. 

@param obj The chirality to be compared to.
@return -# if this is less than obj, 0 if equal, +# if this is greater than obj
*/
int Chirality::compare(Chirality &obj)
{
	return n * 100 + m - (obj.getn() * 100 + obj.getm());
}

/**
Compares Chiralities
*/
bool Chirality::operator==(const Chirality val)
{
	if (n == val.getn() && m == val.getm())
	{
		return true;
	}
	return false;
}
