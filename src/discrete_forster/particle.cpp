#include <iostream>
#include <limits>

#include "./particle.h"

namespace mc
{
  // perform free flight within the simulation domain
  void particle::fly(double dt, const std::vector<scatterer>& s_list) {

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
      
      // arma::vec t = s_list[next].pos();
      if (next<0){
        std::cout << "next is less than zero!!!" << std::endl;
        std::cout << "_scat_ptr->left: " << _scat_ptr->left << std::endl;
        std::cout << "_scat_ptr->right: " << _scat_ptr->right << std::endl;
      }

      arma::vec t = s_list.at(next).pos();
      double dist = arma::norm(_pos - t);

      if (dist / _velocity < dt) {
        // _pos = s_list[next].pos();
        // _scat_ptr = &(s_list[next]);
        _pos = s_list.at(next).pos();
        _scat_ptr = &(s_list.at(next));
        dt -= dist / _velocity;
      } else {
        // arma::vec d = arma::normalise(s_list[next].pos() - _pos);
        arma::vec d = arma::normalise(s_list.at(next).pos() - _pos);
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
