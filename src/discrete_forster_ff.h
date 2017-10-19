#ifndef discrete_forster_ff_h
#define discrete_forster_ff_h

#include <iostream>
#include <array>

#include "utility.h"
#include "discrete_forster_particle.h"

namespace mc
{

class discrete_forster_free_flight
{
public:
	typedef mc::discrete_forster_particle t_particle; // particle type
	// typedef mc::discrete_forster_free_flight t_ff; // free_flight type
	// typedef mc::discrete_forster_region t_region; // region type
	// typedef mc::discrete_forster_scatter t_scatter; // scatter
private:

public:
	discrete_forster_free_flight() {}; // constructor
	void fly(t_particle* p, const mc::t_float &dt) {}; // perform free_flight
	void check_boundary(t_particle* p, const mc::t_float &dt, const std::pair<mc::arr1d, mc::arr1d>& domain); // check for collision to boundaries
}; // end class discrete_forster_free_flight

} // end namespace mc

#endif // discrete_forster_ff_h
