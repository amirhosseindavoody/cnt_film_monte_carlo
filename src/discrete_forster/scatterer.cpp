#include <iostream>
#include <limits>

#include "discrete_forster_particle.h"
#include "scatterer.h"

namespace mc
{
	void scatterer::update_state(discrete_forster_particle* p) // update the final state of the particle
	{
		mc::t_float dice = _max_rate*double(std::rand())/double(RAND_MAX);
		auto it = _neighbors.begin();
		auto last = std::prev(_neighbors.end());

		while ((it->rate <= dice) and (it!=last))
		{
			it++;
		}

		p->set_pos(it->s_ptr->pos());
		p->set_scatterer(it->s_ptr);
	};
} // mc namespace
