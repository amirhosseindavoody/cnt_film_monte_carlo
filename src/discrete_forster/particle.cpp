#include <iostream>
#include <limits>

#include "./particle.h"

namespace mc
{
  // perform free flight within the simulation domain
  void particle::fly(double dt, const std::pair<arma::vec, arma::vec>& domain) {

    std::cout << "new free fly:" << std::endl;

    std::cout << "1" << std::endl;
    scatterer* next_scat = nullptr;

    while (true){
      if (_scat_ptr == nullptr) {
        std::cout << "_scat_ptr is NULL!!!" << std::endl;
      }

      try{

        // determine next scatterer
        if (_heading_right) {
          if (_scat_ptr->right){
            next_scat = _scat_ptr->right;
          } else {
            next_scat = _scat_ptr->left;
          }
        } else {
          if (_scat_ptr->left) {
            next_scat = _scat_ptr->left;
          } else {
            next_scat = _scat_ptr->right;
          }
        }

      } catch (...){
        std::cout << "failed!!!" << std::endl;
        std::exit(0);
      }

      if (next_scat == nullptr) {
        std::cout << "next_scat is nullptr!!!" << std::endl;
      }

      if (next_scat == NULL) {
        std::cout << "next_scat is NULL!!!" << std::endl;
      }

      // determine the movement direction
      if (next_scat==_scat_ptr->right)
        _heading_right = true;
      else
        _heading_right = false;

      std::cout << "next_scat=" << next_scat << std::endl;
      
      double dist;
      std::cout << "2" << std::endl;
      try {
        dist = arma::norm(_pos - next_scat->pos());
      } catch (...) {
        std::cout << "failed here!!!" << std::endl;
        std::exit(0);
      }

      std::cout << dt << "," << dist << std::endl;

      if (dist / _velocity < dt) {
        _pos = next_scat->pos();
        _scat_ptr = next_scat;
        dt -= dist / _velocity;
      } else {
        arma::vec d = arma::normalise(next_scat->pos() - _pos);
        _pos += (_velocity * dt * d);
        return;
      }

    }

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
