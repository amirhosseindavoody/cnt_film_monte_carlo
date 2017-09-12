#include <iostream>
#include <limits>

#include "ff.h"

namespace mc
{

// constructor
free_flight::free_flight(mc::arr1d accel)
{
	// make sure that there is an inverse for acceleration value. othersize put a small value instead of acceleration so that 1/accel is not infinity
	for (auto &elem : accel)
	{
		if (! std::isfinite(1./elem))
		{
			elem = std::numeric_limits<mc::t_float>::epsilon();
		}
	}

	for (int i=0; i<accel.size(); ++i)
	{
		_acceleration[i] = accel[i];
		_inv_acceleration[i] = 1./accel[i];
	}
};

// perform free flight
void free_flight::fly(mc::arr1d &pos, mc::arr1d &velocity, const mc::t_float &eff_mass, const mc::t_float &dt)
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
void free_flight::check_boundary(mc::arr1d &pos, mc::arr1d &velocity, const mc::arr1d& old_pos, const mc::arr1d& old_velocity, const mc::t_float &eff_mass, const mc::t_float &dt, const std::pair<mc::arr1d, mc::arr1d>& domain)
{
	mc::t_float t_collision;
	mc::t_float temp_expr;

	for (int i=0; i<pos.size(); ++i)
	{
		if (pos[i] > domain.second[i])
		{
			temp_expr = std::abs(2*_acceleration[i]*(old_pos[i]-domain.second[i])/std::pow(old_velocity[i],2));
			if (temp_expr < 1.e-3)
			{
				t_collision = (domain.second[i]-old_pos[i])/old_velocity[i];
			}
			else
			{
				t_collision = (-old_velocity[i]+std::sqrt(std::pow(old_velocity[i],2)-2.*_acceleration[i]*(old_pos[i]-domain.second[i])))*(_inv_acceleration[i]);
			}
			pos[i] = domain.second[i] + (-old_velocity[i]-_acceleration[i]*t_collision)*(dt-t_collision) + _acceleration[i]*std::pow(dt-t_collision,2)/2.;
			velocity[i] = ((-old_velocity[i]-_acceleration[i]*t_collision)+_acceleration[i]*(dt-t_collision));
		}
		else if(pos[i] < domain.first[i])
		{
			temp_expr = std::abs(2*_acceleration[i]*(old_pos[i]-domain.first[i])/std::pow(old_velocity[i],2));
			if (temp_expr < 1.e-3)
			{
				t_collision = (domain.first[i]-old_pos[i])/old_velocity[i];
			}
			else
			{
				t_collision = (-old_velocity[i]-std::sqrt(std::pow(old_velocity[i],2)-2.*_acceleration[i]*(old_pos[i]-domain.first[i])))*(_inv_acceleration[i]);
			}
			pos[i] = domain.first[i] + (-old_velocity[i]-_acceleration[i]*t_collision)*(dt-t_collision) + _acceleration[i]*std::pow(dt-t_collision,2)/2.;
			velocity[i] = ((-old_velocity[i]-_acceleration[i]*t_collision)+_acceleration[i]*(dt-t_collision));
		}
	}
};

} // mc namespace
