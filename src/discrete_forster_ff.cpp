#include <iostream>
#include <array>

#include "utility.h"
#include "discrete_forster_ff.h"

namespace mc
{

void discrete_forster_free_flight::check_boundary(discrete_forster_free_flight::t_particle* p, const mc::t_float &dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // check for collision to boundaries
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

} // end namespace mc
