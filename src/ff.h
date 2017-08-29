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
	void fly(mc::arr1d &pos, mc::arr1d &mom, const mc::t_float &eff_mass, const mc::t_float &dt, const mc::arr1d& volume); // perform free_flight
	const mc::arr1d& acceleration(); // get constant reference to acceleration vector
	void check_boundary(mc::arr1d &pos, mc::arr1d &mom, const mc::arr1d& old_pos, const mc::arr1d& old_mom, const mc::t_float &eff_mass, const mc::t_float &dt, const mc::arr1d& volume); // check for collision to boundaries
}; // end class particle

} // end namespace mc

#endif // ff_h