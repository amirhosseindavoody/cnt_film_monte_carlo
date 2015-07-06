/**
CNT.cpp
Stores all relevant information for a carbon nanotube

@author Alex Gabourie
@version 1.00
*/

#include "stdafx.h"
#include "CNT.h"
#define _USE_MATH_DEFINES
#include <math.h>


/**
Reads a CNT file and creates a CNT object with all the information stored
in that file.

@param filePath The path of the file containing the CNT info
@return CNT Object
*/
CNT::CNT(string filePath)
{
	
}

//Destructor implemented automatically


/**
Uses the chirality of the nanotube to calculate the diameter
@param n Hamada n paramter
@param m Hamada m parameter
*/
void CNT::setDiameter(int n, int m)
{
	diameter = A_CC*sqrt(pow(n*1., 2) + pow(m*1., 2) + n*m) / M_PI;
}

/**
Gets the diameter of the CNT

@param void
@return Diameter of the CNT
*/
double CNT::getDiameter()
{
	return diameter;
}

/**
Gets the length of the CNT

@param void
@return Length of the CNT
*/
double CNT::getLength()
{
	return length;
}

/**
Gets the height of each compositional cylinders

@param void
@return The height of each compositional cylinders
*/
double CNT::getCylinderHeight()
{
	return cylinderHeight;
}

/**
Gets the separation between two compositional cylinders

@param void
@return The separation between two compositional cylinders
*/
double CNT::getTubeSeparation()
{
	return tubeSeparation;
}

/**
Gets the minimum spacing between two nanotubes

@param void
@return minimum spacing between two nanotubes
*/
double CNT::getMinSpacing()
{
	return minSpacing;
}

/**
Gets m parameter of the CNT

@param void
@return m
*/
int CNT::getm()
{
	return m;
}

/**
Gets n parameter of the CNT

@param void
@return n
*/
int CNT::getn()
{
	return n;
}




