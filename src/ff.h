#ifndef ff_h
#define ff_h

#include <iostream>
#include <array>

#include "utility.h"

namespace mc
{

class free_flight
{
private:
	mc::arr1d _acceleration; // acceleration in the three dimension space

public:
	free_flight(mc::arr1d accel = {0,0,0}); // constructor
	void fly(mc::arr1d &pos, mc::arr1d &mom, const mc::t_float &eff_mass, const mc::t_float &dt); // perform free_flight

}; // end class particle

} // end namespace mc

#endif // ff_h