#include <iostream>
#include <limits>

#include "discrete_forster_particle.h"

namespace mc
{
  //constructor
	discrete_forster_particle::discrete_forster_particle(const mc::arr1d& pos, const std::shared_ptr<t_ff>& pilot, const std::shared_ptr<t_scatter>& m_scatterer) // constructor
	{
		set_pos(pos);
		set_old_pos(pos); // we set the old position separately so that it is not filled with some non-sense
		set_pilot(pilot);
		set_scatterer(m_scatterer);
		_ff_time = scatterer()->ff_time();
	};

  // perform free flight within the simulation domain
  void discrete_forster_particle::fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // perform free flight within the simulation domain
  {
    pilot()->check_boundary(this, dt, domain);
  };

  // reinitialize particle properties instead of creating new particles
	void discrete_forster_particle::reinitialize(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner, const mc::t_float& beta, const mc::t_float& mass, const std::shared_ptr<t_ff>& pilot, const std::shared_ptr<t_scatter>& m_scatterer)
	{
		for (int i=0; i<lower_corner.size(); ++i)
		{
			set_pos(i, lower_corner[i]+(upper_corner[i]-lower_corner[i])*mc::get_rand_include_zero<mc::t_float>());
		}
		_ff_time = scatterer()->ff_time();
	};

  void discrete_forster_particle::step(mc::t_float dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // step particle state for dt in time
	{
		while(_ff_time <= dt)
		{
			dt -= _ff_time;
			fly(_ff_time, domain);
			_scatterer->update_state(this);
			_ff_time = _scatterer->ff_time();
		}
		fly(dt, domain);
		_ff_time -= dt;
	};

  // update the _ff_time by calling the underlying scatterer
	void discrete_forster_particle::get_ff_time()
  {
    _ff_time = _scatterer->ff_time();
  };

} // mc namespace
