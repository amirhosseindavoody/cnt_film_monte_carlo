#pragma once

#include <Eigen>

#include "CNT.h"
#include "segment.h"

using namespace std;
using namespace Eigen;

class tableElem
{ 

	double r; //distance between segment and segment with ID: tubeidx, segidx
	double theta; //Angle in radians between two segments
	double gamma; //distance^6/time constant for initial purposes
	double gammaTot; //Transistion rate to be used in simulation
	int tubeidx; //tube interacting with from current tube
	int segidx; //segment number on the tube of tubeidx

private:
	void setRate();
	void setr(double rnew);
	void setTheta(double t);
	void setGamma(double g);
	void setTubeidx(int num);
	void setSegidx(int num);

public:
	tableElem(double rnew=1, double t=0, double g=0, int tube=0, int seg=0);
	double getRate();
	double getr();
	double getTheta();
	double getGamma();
	int getTubeidx();
	int getSegidx();
	static double calcDist(Vector3d v1, Vector3d v2);
	static double calcThet(segment &s1, segment &s2);
};

