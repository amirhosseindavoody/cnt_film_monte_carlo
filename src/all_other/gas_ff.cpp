#include <iostream>
#include <limits>

#include "particle.h"
#include "gas_ff.h"

namespace mc
{
  // perform free_flight
	void gas_free_flight::fly(mc::particle* p, const mc::t_float& dt)
	{
		const mc::t_float coeff = dt*dt/2.;
		for (int i=0; i<p->pos().size(); ++i)
		{
			p->set_pos(i, p->pos(i)+p->velocity(i)*dt + _acceleration[i]*coeff);
			p->set_velocity(i, p->velocity(i)+_acceleration[i]*dt);
		}
	};

} // mc namespace
