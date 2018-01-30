#ifndef forster_scatter_h
#define forster_scatter_h

#include <iostream>
#include <array>

#include "../helper/utility.h"
#include "scatter.h"

namespace mc
{

class forster_scatter : public scatter
{
private:
	mc::t_float _max_rate; // maximum scattering rate in the scattering table
	mc::t_float _inverse_max_rate; // inverse of the maximum scattering rate which is the lifetime
	mc::t_float _min_radius; // minimum distance between hopping sites in the forster scattering process

public:
	inline forster_scatter() // constructor
	{
		_max_rate = 1.e14;
		_inverse_max_rate = 1./_max_rate;
		_min_radius = 5.e-10; // in meters units
	};
	inline mc::t_float ff_time() const // get random free flight time
	{
		return -_inverse_max_rate*std::log(mc::get_rand_exclude_zero<mc::t_float>());
	};
	inline void update_state(mc::particle* p) // update the final state of the particle
	{
		// get uniformly distribution direction
		mc::t_float theta = std::acos(1.-2.*mc::get_rand_include_zero<mc::t_float>());
		mc::t_float phi = 2.*mc::pi*mc::get_rand_include_zero<mc::t_float>();
		mc::t_float radius = _min_radius/std::cbrt(mc::get_rand_exclude_zero<mc::t_float>());
		mc::arr1d new_pos = {p->pos(0)+radius*std::sin(theta)*std::cos(phi), p->pos(1)+radius*std::sin(theta)*std::sin(phi), p->pos(2)+radius*std::cos(theta)};
		p->set_pos(new_pos);
	};

}; // end class scatter

} // end namespace mc

#endif // gas_scatter_h
