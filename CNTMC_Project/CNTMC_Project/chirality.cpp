#include "chirality.h"


/**
Sets the Chirality object to some default values. DO NOT USE CNTs CONSTRUCTED
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
int Chirality::getm()
{
	return m;
}

/**
Gets n parameter of the CNT

@param void
@return n
*/
int Chirality::getn()
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
	m = m_new;
}