#include <iostream>
#include <limits>

#include "discrete_forster_scatter.h"
#include "discrete_forster_particle.h"

namespace mc
{
	void discrete_forster_scatter::update_state(discrete_forster_scatter::t_particle* p) // update the final state of the particle
	{
		mc::t_float dice = _max_rate*mc::get_rand_include_zero<mc::t_float>();
		auto it = _neighbors.begin();
		auto last = std::prev(_neighbors.end());
		while ((std::get<2>(*it) <= dice) and (it!=last))
		{
			it++;
		}

		p->set_pos(std::get<0>(*it)->pos());
		p->set_scatterer(std::get<0>(*it));
	};
} // mc namespace
