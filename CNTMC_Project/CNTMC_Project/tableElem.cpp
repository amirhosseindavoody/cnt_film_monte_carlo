#include "stdafx.h"
#include "tableElem.h"


/**
Creates table element object

@return tableElem Object
*/
tableElem::tableElem()
{
	
}

/**
Destructor for class
*/
tableElem::~tableElem()
{
}

/**
Sets the distance r

@param rnew Distance between one seg to other
*/
void tableElem::setr(double rnew)
{
	if (rnew < 0)
	{
		cout << "Error: Negative r not accepted.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	r = rnew;
}

/**
Sets the angle between two segments

@param t The angle between two segments
*/
void tableElem::setTheta(double t)
{
	if (t < 0 || t > 2*M_PI)
	{
		cout << "Error: Theta must be between 0 and 2*pi.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	theta = t;
}

/**
Sets gamma parameter

@param g The new gamma value
*/
void tableElem::setGamma(double g)
{
	if (g < 0)
	{
		cout << "Error: g in setGamma must be positive.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	gamma = g;
}

/**
Sets the destination tube number

@param num The tube number
*/
void tableElem::setTubeNum(int num)
{
	if (num < 0)
	{
		cout << "Error: Tube number must not be negative.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	tubeNum = num;
}

/**
Sets the destination segment number

@param num The segment number
*/
void tableElem::setSegNum(int num)
{
	if (num < 0)
	{
		cout << "Error: Segment number must not be negative.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	segNum = num;
}

/**
Gets the total transition rate based on gamma, r, and theta

@return The transition rate in inverse seconds
*/
double tableElem::getRate()
{
	return gamma*cos(theta) / pow(r, 6);
}

/**
Sets the distance r

@param rnew Distance between one seg to other
*/
double tableElem::calcDist(Vector3d v1, Vector3d v2)
{
	return (v1 - v2).norm();
}

/**
Sets the distance r

@param rnew Distance between one seg to other
*/
double tableElem::calcThet(Vector3d v1, Vector3d v2)
{
	return acos(v1.dot(v2) / (v1.norm()*v2.norm()));
}


