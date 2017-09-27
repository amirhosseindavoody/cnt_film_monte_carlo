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
	void fly(mc::particle* p, const mc::t_float &dt) {}; // perform free_flight
	void check_boundary(mc::particle* p, const mc::t_float &dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // check for collision to boundaries
	{
		for (int i=0; i<domain.first.size(); ++i)
		{
			if (p->pos(i) < domain.first[i])
			{
				p->rewind_pos();
				return;
			}
			if (p->pos(i) > domain.second[i])
			{
				p->rewind_pos();
				return;
			}
		}
	};
}; // end class particle

} // end namespace mc

#endif // forster_ff_h
