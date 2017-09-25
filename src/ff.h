#ifndef ff_h
#define ff_h

#include "utility.h"

namespace mc
{
class particle;

class free_flight
{
private:

public:
	virtual void fly(mc::particle* p, const mc::t_float& dt) = 0; // perform free_flight on particle
	virtual void check_boundary(mc::particle* p, const mc::t_float &dt, const std::pair<mc::arr1d, mc::arr1d>& domain) = 0; // check for collision to boundaries
}; // end class particle

} // end namespace mc

#endif // ff_abstract_h
