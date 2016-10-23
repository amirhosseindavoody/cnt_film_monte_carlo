#include <stdio.h>
#include <iostream>

#include "tableElem.h"

// Creates table element object
tableElem::tableElem(double rnew, double t, double g, int tube, int seg)
{
	r = rnew;
	theta = t;
	gamma = g;
	tubeidx = tube;
	segidx = seg;

	setRate();
}

void tableElem::setRate()
{
	gammaTot = abs(gamma*cos(theta) / pow(r, 6));
}


// Sets the distance r
void tableElem::setr(double rnew)
{
	if (rnew < 0)
	{
		cout << "Error: Negative r not accepted." << endl;
		exit(EXIT_FAILURE);
	}
	r = rnew;
}


// Sets the angle between two segments
void tableElem::setTheta(double t)
{
	if (t < 0 || t > 2*M_PI)
	{
		cout << "Error: Theta must be between 0 and 2*pi." << endl;
		exit(EXIT_FAILURE);
	}
	theta = t;
}


// Sets gamma parameter
void tableElem::setGamma(double g)
{
	if (g < 0)
	{
		cout << "Error: g in setGamma must be positive." << endl;
		exit(EXIT_FAILURE);
	}
	gamma = g;
}


// Sets the destination tube number
void tableElem::setTubeidx(int num)
{
	if (num < 0)
	{
		cout << "Error: Tube number must not be negative." << endl;
		exit(EXIT_FAILURE);
	}
	tubeidx = num;
}


// Sets the destination segment number
void tableElem::setSegidx(int num)
{
	if (num < 0)
	{
		cout << "Error: Segment number must not be negative." << endl;
		exit(EXIT_FAILURE);
	}
	segidx = num;
}


// Gets the total transition rate based on gamma, r, and theta
double tableElem::getRate()
{
	return gammaTot;
}


// Gets r value
double tableElem::getr()
{
	return r;
}


// Gets theta value
double tableElem::getTheta()
{
	return theta;
}


// Gets gamma value
double tableElem::getGamma()
{
	return gamma;
}

// Gets tube number
int tableElem::getTubeidx()
{
	return tubeidx;
}


// Gets the segment number
int tableElem::getSegidx()
{
	return segidx;
}



// Calculates distance between two segments
double tableElem::calcDist(Vector3d v1, Vector3d v2)
{
	return (v1 - v2).norm();
}


// Calculates the angle between two vectors
double tableElem::calcThet(segment &s1, segment &s2)
{
	Vector3d v1 = s1.p2 - s1.p1;
	Vector3d v2 = s2.p2 - s2.p1;
	double val = acos(v1.dot(v2) / (v1.norm()*v2.norm())); //range 0 to pi
	if (val <= M_PI / 2.0)
	{
		return val;
	}
	return (M_PI - val);

}


