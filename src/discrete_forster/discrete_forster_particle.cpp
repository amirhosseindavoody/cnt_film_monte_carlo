#include <iostream>
#include <limits>

#include "discrete_forster_particle.h"

namespace mc
{
  //constructor
  discrete_forster_particle::discrete_forster_particle(
      const arma::vec& pos, const std::shared_ptr<t_ff>& pilot,
      const std::shared_ptr<scatterer>& s) {
  set_pos(pos);
  // initialize old position with the current position just for safety
  set_old_pos(pos); 
  set_pilot(pilot);
  set_scatterer(s);
  _ff_time = scat_ptr()->ff_time();
	};

  // perform free flight within the simulation domain
  void discrete_forster_particle::fly(
      const double& dt, const std::pair<arma::vec, arma::vec>& domain) {
    pilot()->check_boundary(this, dt, domain);
  };

  // reinitialize particle properties instead of creating new particles
  void discrete_forster_particle::reinitialize(
      const arma::vec& lower_corner, const arma::vec& upper_corner,
      const double& beta, const double& mass,
      const std::shared_ptr<t_ff>& pilot,
      const std::shared_ptr<scatterer>& m_scatterer) {
    set_pos(lower_corner +
            ((upper_corner - lower_corner) % arma::randu<arma::vec>(3)));
    _ff_time = scat_ptr()->ff_time();
	};

  // step particle state for dt in time
  void discrete_forster_particle::step(
      double dt, const std::pair<arma::vec, arma::vec>& domain) {
    while (_ff_time <= dt) {
      dt -= _ff_time;
      fly(_ff_time, domain);
      _scat_ptr->update_state(this);
      _ff_time = _scat_ptr->ff_time();
    }
    fly(dt, domain);
    _ff_time -= dt;
	};

  // update the _ff_time by calling the underlying scatterer
	void discrete_forster_particle::get_ff_time() {
    _ff_time = _scat_ptr->ff_time();
  };

} // mc namespace
