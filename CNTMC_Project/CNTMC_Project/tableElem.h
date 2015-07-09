#pragma once

#include <Eigen>
#include "CNT.h"

using namespace std;
using namespace Eigen;

struct segment;

class tableElem
{ 

	double r; //distance between segment and segment with ID: tubeNum, segNum
	double theta; //Angle in radians between two segments
	double gamma; //distance^6/time constant for initial purposes
	double gammaTot; //Transistion rate to be used in simulation
	int tubeNum; //tube interacting with from current tube
	int segNum; //segment number on the tube of tubeNum

private:
	void setRate();
	void setr(double rnew);
	void setTheta(double t);
	void setGamma(double g);
	void setTubeNum(int num);
	void setSegNum(int num);

public:
	tableElem();
	tableElem(double rnew, double t, double g, int tube, int seg);
	~tableElem();
	double getRate();
	double getr();
	double getTheta();
	double getGamma();
	int getTubeNum();
	int getSegNum();
	static double calcDist(Vector3d v1, Vector3d v2);
	static double calcThet(vector<segment>::iterator s1, vector<segment>::iterator s2);
};

