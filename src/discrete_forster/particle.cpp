#include <iostream>
#include <limits>

#include "./particle.h"

namespace mc
{
  // perform free flight within the simulation domain
  void particle::fly(double dt, const std::vector<scatterer>& s_list) {

    if (_scat_ptr->left<0 && _scat_ptr->right<0)
      return;

    int next;

    while (true){

      // determine next scatterer
      if (_heading_right) {
        if (_scat_ptr->right>-1){
          next = _scat_ptr->right;
        } else {
          next = _scat_ptr->left;
        }
      } else {
        if (_scat_ptr->left>-1) {
          next = _scat_ptr->left;
        } else {
          next = _scat_ptr->right;
        }
      }

      // determine the movement direction
      if (next==_scat_ptr->right)
        _heading_right = true;
      else
        _heading_right = false;

      std::valarray<double> t = s_list[next].pos();
      double dist = std::sqrt(std::pow(_pos - t, 2).sum());

      if (dist / _velocity < dt) {
        _pos = s_list[next].pos();
        _scat_ptr = &(s_list[next]);
        dt -= dist / _velocity;
      } else {
        std::valarray<double> d = (s_list[next].pos() - _pos);
        d /= std::sqrt((d*d).sum());
        _pos += (_velocity * dt * d);
        return;
      }

    }

  }

  // step particle state for dt in time
  void particle::step(double dt, const std::vector<scatterer>& s_list, const double& max_hop_radius) {
    
    _old_pos = _pos;
    const scatterer* new_scat_ptr = nullptr;

    while (_ff_time <= dt) {
      dt -= _ff_time;
      
      fly(_ff_time, s_list);
      
      new_scat_ptr = _scat_ptr->update_state(max_hop_radius);
      
      if (new_scat_ptr!=_scat_ptr){
        _scat_ptr = new_scat_ptr;
        _pos = _scat_ptr->pos();
      }

      _ff_time = _scat_ptr->ff_time();
    }
    
    fly(dt, s_list);

    _ff_time -= dt;
	};

} // mc namespace
