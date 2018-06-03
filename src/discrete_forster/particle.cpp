#include <iostream>
#include <limits>

#include "./particle.h"

namespace mc
{
  // perform free flight within the simulation domain
  void particle::fly( const double& dt, const std::pair<arma::vec, arma::vec>& domain) {

    arma::vec prev_pos = _pos;

    // do some free fly stuff

    // rewind the position if we are out of boundary
    if (arma::any(_pos < domain.first) || arma::any(_pos > domain.second)) {
      _pos = prev_pos;
    }
  };

  // step particle state for dt in time
  void particle::step(double dt, const std::pair<arma::vec, arma::vec>& domain, const double& max_hop_radius) {
    _old_pos = _pos;
    while (_ff_time <= dt) {
      dt -= _ff_time;
      // fly(_ff_time, domain);
      _scat_ptr = _scat_ptr->update_state(max_hop_radius);
      _pos = _scat_ptr->pos();
      _ff_time = _scat_ptr->ff_time();
    }
    // fly(dt, domain);
    _ff_time -= dt;
	};

} // mc namespace
