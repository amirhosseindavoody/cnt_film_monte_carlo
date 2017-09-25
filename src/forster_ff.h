#ifndef forster_ff_h
#define forster_ff_h

#include <iostream>
#include <array>

#include "utility.h"

namespace mc
{

class forster_free_flight : public free_flight
{
private:

public:
	forster_free_flight() {}; // constructor
	void fly(mc::t_float* p, const mc::t_float &dt) {}; // perform free_flight
	void check_boundary(mc::particle* p, const mc::t_float &dt, const std::pair<mc::arr1d, mc::arr1d>& domain) {}; // check for collision to boundaries
}; // end class particle

} // end namespace mc

#endif // forster_ff_h
