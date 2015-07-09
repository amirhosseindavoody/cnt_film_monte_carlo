#include "stdafx.h"
#include "tableElem.h"
#include <iostream>


tableElem::tableElem()
{
	setr(1);
	setTheta(0);
	setGamma(0);
	setRate();
	setTubeNum(0);
	setSegNum(0);
}
/**
Creates table element object

@return tableElem Object
*/
tableElem::tableElem(double rnew, double t, double g, int tube, int seg)
{
	setr(rnew);
	setTheta(t);
	setGamma(g);
	setRate();
	setTubeNum(tube);
	setSegNum(seg);
}

/**
Destructor for class
*/
tableElem::~tableElem()
{
}

void tableElem::setRate()
{
	gammaTot = gamma*cos(theta) / pow(r, 6);
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
	return gammaTot;
}

/**
Gets r value

@return r value
*/
double tableElem::getr()
{
	return r;
}

/**
Gets theta value

@return theta value
*/
double tableElem::getTheta()
{
	return theta;
}

/**
Gets gamma value

@return gamma value
*/
double tableElem::getGamma()
{
	return gamma;
}

/**
Gets tube number

@return tube number
*/
int tableElem::getTubeNum()
{
	return tubeNum;
}

/**
Gets the segment number

@return segment number
*/
int tableElem::getSegNum()
{
	return segNum;
}


/**
Calculates distance between two segments

@param v1 first segment center
@param v2 second segment center

@return The distance between v1 and v2
*/
double tableElem::calcDist(Vector3d v1, Vector3d v2)
{
	return (v1 - v2).norm();
}

/**
Calculates the angle between two vectors

@param s1 first segment
@param s2 second segment
*/
double tableElem::calcThet(vector<segment>::iterator s1, vector<segment>::iterator s2)
{
	Vector3d v1 = s1->p2 - s1->p1;
	Vector3d v2 = s2->p2 - s2->p1;
	return acos(v1.dot(v2) / (v1.norm()*v2.norm())); //range 0 to pi

}


