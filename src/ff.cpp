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
void free_flight::fly(mc::arr1d &pos, mc::arr1d &velocity, const mc::t_float &eff_mass, const mc::t_float &dt, const mc::arr1d& volume)
{
	const mc::t_float coeff = dt*dt/2.;

	for (int i=0; i<pos.size(); ++i)
	{
		pos[i] +=  velocity[i]*dt + _acceleration[i]*coeff;
		velocity[i] += _acceleration[i]*dt;
	}
};

// get constant reference to acceleration vector
const mc::arr1d& free_flight::acceleration()
{
	return _acceleration;
};

// check for collision to boundaries
void free_flight::check_boundary(mc::arr1d &pos, mc::arr1d &velocity, const mc::arr1d& old_pos, const mc::arr1d& old_velocity, const mc::t_float &eff_mass, const mc::t_float &dt, const mc::arr1d& volume)
{
	mc::t_float t_collision;
	const mc::t_float coeff = dt*dt/2.;

	for (int i=0; i<pos.size(); ++i)
	{
		if (pos[i] > volume[i])
		{
			t_collision = (-old_velocity[i]+std::sqrt(std::pow(old_velocity[i],2)-2.*_acceleration[i]*(old_pos[i]-volume[i])))/(_acceleration[i]);
			pos[i] = volume[i] + (-old_velocity[i]-_acceleration[i]*t_collision)*(dt-t_collision) + _acceleration[i]*std::pow(dt-t_collision,2)/2.;
			velocity[i] = ((-old_velocity[i]-_acceleration[i]*t_collision)+_acceleration[i]*(dt-t_collision));

			std::cout << "t_collision = " << t_collision << std::endl;
			std::cin.ignore();
		}
		else if(pos[i] < 0)
		{
			t_collision = (-old_velocity[i]-std::sqrt(std::pow(old_velocity[i],2)-2.*_acceleration[i]*(old_pos[i])))/(_acceleration[i]);
			pos[i] = (-old_velocity[i]-_acceleration[i]*t_collision)*(dt-t_collision) + _acceleration[i]*std::pow(dt-t_collision,2)/2.;
			velocity[i] = ((-old_velocity[i]-_acceleration[i]*t_collision)+_acceleration[i]*(dt-t_collision));

			std::cout << "t_collision = " << t_collision << std::endl;
			std::cin.ignore();
		}
	}
}; 

} // mc namespace