#include <iostream>
#include <limits>

#include "./particle.h"

namespace mc
{
  // perform free flight within the simulation domain
  void particle::fly( const double& dt, const std::pair<arma::vec, arma::vec>& domain) {

    // arma::vec prev_pos = _pos;

    scatterer* next_scat = nullptr;

    if (_heading_right) {
      if (_scat_ptr->right){
        next_scat = _scat_ptr->right;
      } else {
        next_scat = _scat_ptr->left;
        _heading_right = false;
      }
    } else {
      if (_scat_ptr->left) {
        next_scat = _scat_ptr->left;
      } else {
        next_scat = _scat_ptr->right;
        _heading_right = true;
      }
    }

    arma::vec d = arma::normalise(next_scat->pos() - _pos);
    _pos += (_velocity * dt * d);
    if (arma::dot(_scat_ptr->pos() - _pos, next_scat->pos() - _pos) >= 0) {
      _scat_ptr = next_scat;
    }

    // // rewind the position if we are out of boundary
    // if (arma::any(_pos < domain.first) || arma::any(_pos > domain.second)) {
    //   _pos = prev_pos;
    // }
  }

  // step particle state for dt in time
  void particle::step(double dt, const std::pair<arma::vec, arma::vec>& domain, const double& max_hop_radius) {
    
    _old_pos = _pos;
    const scatterer* new_scat_ptr = nullptr;

    while (_ff_time <= dt) {
      dt -= _ff_time;
      
      fly(_ff_time, domain);
      
      new_scat_ptr = _scat_ptr->update_state(max_hop_radius);
      
      if (new_scat_ptr!=_scat_ptr){
        _scat_ptr = new_scat_ptr;
        _pos = _scat_ptr->pos();
      }

      _ff_time = _scat_ptr->ff_time();
    }
    
    fly(dt, domain);

    _ff_time -= dt;
	};

} // mc namespace
