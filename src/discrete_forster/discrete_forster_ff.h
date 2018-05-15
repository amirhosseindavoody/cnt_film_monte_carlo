#ifndef discrete_forster_ff_h
#define discrete_forster_ff_h

#include <iostream>
#include <array>

#include "../helper/utility.h"
#include "./discrete_forster_particle.h"

namespace mc
{

class discrete_forster_free_flight
{
public:
	typedef mc::discrete_forster_particle t_particle; // particle type
private:

public:
	discrete_forster_free_flight() {}; // constructor
	void fly(t_particle* p, const double &dt) {}; // perform free_flight
	void check_boundary(t_particle* p, const double &dt, const std::pair<arma::vec, arma::vec>& domain); // check for collision to boundaries
}; // end class discrete_forster_free_flight

} // end namespace mc

#endif // discrete_forster_ff_h
