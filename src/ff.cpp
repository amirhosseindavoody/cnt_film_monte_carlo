#include <iostream>

#include "ff.h"

namespace mc
{

// constructor
free_flight::free_flight(mc::arr1d accel)
{
	_acceleration = accel;
};

// perform free flight
void free_flight::fly(mc::arr1d &pos, mc::arr1d &mom, const mc::t_float &eff_mass, const mc::t_float &dt)
{
	for (int i=0; i<pos.size(); ++i)
	{
		pos[i] = pos[i] + mom[i]*dt/eff_mass + _acceleration[i]*(dt*dt)/2.;
		mom[i] = mom[i] + _acceleration[i]*dt;
	}
};

} // mc namespace