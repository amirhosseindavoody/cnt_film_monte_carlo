#include "stdafx.h"
#include "tableElem.h"
#include <iostream>

/**
Creates table element object

@return tableElem Object
*/
tableElem::tableElem()
{
	setr(1);
	setTheta(0);
	setGamma(0);
	setRate();
	setTubeidx(0);
	setSegidx(0);
}
/**
Creates table element object for building table through calculation

@return tableElem Object
*/
tableElem::tableElem(double rnew, double t, double g, int tube, int seg)
{
	setr(rnew);
	setTheta(t);
	setGamma(g);
	setRate();
	setTubeidx(tube);
	setSegidx(seg);
}

/**
Creates table element for building table through transfer rate table

@return tableElem Object
*/
tableElem::tableElem(double rnew, double t, int tube, int seg, energy new_e, double rate)
{
	setr(rnew);
	setTheta(t);
	setTubeidx(tube);
	setSegidx(seg);
	setEnergy(new_e);
	setRate(rate);
}


/**
Destructor for class
*/
tableElem::~tableElem()
{
}

/**
Sets the energy of the resulting state

@param e_new The energy to set the state at
*/
void tableElem::setEnergy(energy e_new)
{
	e = e_new;
}



/**
Calculates and sets the rate of transfer for exciton following
transfer rates of gamma/r^6 * cos(theta)
*/
void tableElem::setRate()
{
	gammaTot = abs(gamma*cos(theta) / pow(r, 6));
}

/**
Set the rate of transfer for exciton following rate based on 
Amirhossein's calculations
*/
void tableElem::setRate(double new_rate)
{
	if (new_rate <= 0)
	{
		cout << "Error: Must set rate to be positive number.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	gammaTot = new_rate;
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
void tableElem::setTubeidx(int num)
{
	if (num < 0)
	{
		cout << "Error: Tube number must not be negative.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	tubeidx = num;
}

/**
Sets the destination segment number

@param num The segment number
*/
void tableElem::setSegidx(int num)
{
	if (num < 0)
	{
		cout << "Error: Segment number must not be negative.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	segidx = num;
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
int tableElem::getTubeidx()
{
	return tubeidx;
}

/**
Gets the segment number

@return segment number
*/
int tableElem::getSegidx()
{
	return segidx;
}

/**
Gets the energy of the state the table element represents
*/
energy tableElem::getEnergy()
{
	return e;
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
@return The angle between the two segments (0 to pi/2)
*/
double tableElem::calcThet(vector<shared_ptr<segment>>::iterator s1, vector<shared_ptr<segment>>::iterator s2)
{
	Vector3d v1 = (*s1)->p2 - (*s1)->p1;
	Vector3d v2 = (*s2)->p2 - (*s2)->p1;
	double val = acos(v1.dot(v2) / (v1.norm()*v2.norm())); //range 0 to pi
	if (val <= M_PI / 2.0)
	{
		return val;
	}
	return (M_PI - val);

}


