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
void free_flight::fly(mc::arr1d &pos, mc::arr1d &mom, const mc::t_float &eff_mass, const mc::t_float &dt, const mc::arr1d& volume)
{
	const mc::t_float coeff_1 = dt/eff_mass;
	const mc::t_float coeff_2 = dt*dt/2.;

	for (int i=0; i<pos.size(); ++i)
	{
		pos[i] = pos[i] + mom[i]*coeff_1 + _acceleration[i]*coeff_2;
		mom[i] = mom[i] + _acceleration[i]*dt;
	}
};

// get constant reference to acceleration vector
const mc::arr1d& free_flight::acceleration()
{
	return _acceleration;
};

// check for collision to boundaries
void free_flight::check_boundary(mc::arr1d &pos, mc::arr1d &mom, const mc::arr1d& old_pos, const mc::arr1d& old_mom, const mc::t_float &eff_mass, const mc::t_float &dt, const mc::arr1d& volume)
{
	mc::t_float t0;
	const mc::t_float coeff_1 = dt/eff_mass;
	const mc::t_float coeff_2 = dt*dt/2.;

	for (int i=0; i<pos.size(); ++i)
	{
		if (pos[i] > volume[i])
		{
			t0 = 
		}
		else if(pos[i] < 0)
		{

		}
	}
}; 

} // mc namespace