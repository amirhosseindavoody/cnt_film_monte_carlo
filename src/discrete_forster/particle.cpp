#include <iostream>
#include <limits>

#include "./particle.h"

namespace mc
{
  // perform free flight within the simulation domain
  void particle::fly( const double& dt, const std::pair<arma::vec, arma::vec>& domain) {

    if (arma::any(_pos < domain.first)) {
      rewind_pos();
    }

    if (arma::any(_pos > domain.second)) {
      rewind_pos();
    }
  };

  // reinitialize particle properties instead of creating new particles
  void particle::reinitialize(
      const arma::vec& lower_corner, const arma::vec& upper_corner,
      const double& beta, const double& mass,
      const std::shared_ptr<free_flight>& pilot,
      const std::shared_ptr<scatterer>& m_scatterer) {
    set_pos(lower_corner +
            ((upper_corner - lower_corner) % arma::randu<arma::vec>(3)));
    _ff_time = scat_ptr()->ff_time();
	};

  // step particle state for dt in time
  void particle::step(
      double dt, const std::pair<arma::vec, arma::vec>& domain, const double& max_hop_radius) {
    while (_ff_time <= dt) {
      dt -= _ff_time;
      fly(_ff_time, domain);
      _scat_ptr = _scat_ptr->update_state(max_hop_radius);
      _pos = _scat_ptr->pos();
      _ff_time = _scat_ptr->ff_time();
    }
    fly(dt, domain);
    _ff_time -= dt;
	};

} // mc namespace
