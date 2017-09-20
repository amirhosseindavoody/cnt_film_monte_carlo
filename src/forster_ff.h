#ifndef forster_ff_h
#define forster_ff_h

#include <iostream>
#include <array>

#include "utility.h"

namespace mc
{

class forster_free_flight : 
{
private:
	mc::arr1d _acceleration; // acceleration in the three dimension space
	mc::arr1d _inv_acceleration; // inverse of the acceleration in three dimension space

public:
	inline free_flight(mc::arr1d accel = {0,0,0}) // constructor
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
	inline void fly(mc::arr1d &pos, mc::arr1d &velocity, const mc::t_float &eff_mass, const mc::t_float &dt) // perform free_flight
	{
		const mc::t_float coeff = dt*dt/2.;

		for (int i=0; i<pos.size(); ++i)
		{
			pos[i] +=  velocity[i]*dt + _acceleration[i]*coeff;
			velocity[i] += _acceleration[i]*dt;
		}
	};
	inline const mc::arr1d& acceleration() // get constant reference to acceleration vector
	{
		return _acceleration;
	};
	inline void check_boundary(mc::arr1d &pos, mc::arr1d &velocity, const mc::arr1d& old_pos, const mc::arr1d& old_velocity, const mc::t_float &eff_mass, const mc::t_float &dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // check for collision to boundaries
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
}; // end class particle

} // end namespace mc

#endif // forster_ff_h
