#ifndef energy_dispersion_h
#define energy_dispersion_h

#include <iostream>
#include <vector>
#include <memory>

#include "utility.h"

namespace mc
{

class energy_dispersion
{
private:
	mc::v3d _energy; // energy dispersion curves
	mc::v3d _k1, _k2, _k3; // k_vectors over the grid points
	mc::t_int _n1, _n2, _n3; // number of elements in x, y, and z direction

public:
	energy_dispersion(); // constructor

}; //energy_dispersion class

} //mc namespace

#endif // energy_dispersion_h