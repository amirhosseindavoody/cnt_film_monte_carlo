#ifndef ff_gas_h
#define ff_gas_h

#include <iostream>
#include <array>

#include "particle.h"
#include "ff.h"
#include "../helper/utility.h"

namespace mc
{

// class particle;

class gas_free_flight: public free_flight
{
private:
	mc::arr1d _acceleration; // acceleration in the three dimension space
	mc::arr1d _inv_acceleration; // inverse of the acceleration in three dimension space

public:
	gas_free_flight(mc::arr1d accel = {0,0,0}) // constructor
	{
		// make sure that there is an inverse for acceleration value. othersize put a small value instead of acceleration so that 1/accel is not infinity
		for (auto &elem : accel)
		{
			if (! std::isfinite(1./elem))
			{
				elem = std::numeric_limits<mc::t_float>::epsilon();
			}
		}

		for (unsigned i=0; i<accel.size(); ++i)
		{
			_acceleration[i] = accel[i];
			_inv_acceleration[i] = 1./accel[i];
		}
	};
	void fly(mc::particle* p, const mc::t_float& dt); // perform free_flight
	inline const mc::arr1d& acceleration() // get constant reference to acceleration vector
	{
		return _acceleration;
	};
	inline void check_boundary(mc::particle* p, const mc::t_float &dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // check for collision to boundaries
	{
		mc::t_float t_collision;
		mc::t_float temp_expr;

		for (unsigned i=0; i<p->pos().size(); ++i)
		{
			if (p->pos(i) > domain.second[i])
			{
				temp_expr = std::abs(2*_acceleration[i]*(p->old_pos(i)-domain.second[i])/std::pow(p->old_velocity(i),2));
				if (temp_expr < 1.e-3)
				{
					t_collision = (domain.second[i]-p->old_pos(i))/p->old_velocity(i);
				}
				else
				{
					t_collision = (-p->old_velocity(i)+std::sqrt(std::pow(p->old_velocity(i),2)-2.*_acceleration[i]*(p->old_pos(i)-domain.second[i])))*(_inv_acceleration[i]);
				}
				p->set_pos(i, domain.second[i] + (-p->old_velocity(i)-_acceleration[i]*t_collision)*(dt-t_collision) + _acceleration[i]*std::pow(dt-t_collision,2)/2.);
				p->set_velocity(i, (-p->old_velocity(i)-_acceleration[i]*t_collision)+_acceleration[i]*(dt-t_collision));
			}
			else if(p->pos(i) < domain.first[i])
			{
				temp_expr = std::abs(2*_acceleration[i]*(p->old_pos(i)-domain.first[i])/std::pow(p->old_velocity(i),2));
				if (temp_expr < 1.e-3)
				{
					t_collision = (domain.first[i]-p->old_pos(i))/p->old_velocity(i);
				}
				else
				{
					t_collision = (-p->old_velocity(i)-std::sqrt(std::pow(p->old_velocity(i),2)-2.*_acceleration[i]*(p->old_pos(i)-domain.first[i])))*(_inv_acceleration[i]);
				}
				p->set_pos(i, domain.first[i] + (-p->old_velocity(i)-_acceleration[i]*t_collision)*(dt-t_collision) + _acceleration[i]*std::pow(dt-t_collision,2)/2.);
				p->set_velocity(i, (-p->old_velocity(i)-_acceleration[i]*t_collision)+_acceleration[i]*(dt-t_collision));
			}
		}
	};
}; // end class gas_free_flight

} // end namespace mc

#endif // ff_gas_h
