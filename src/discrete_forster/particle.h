#ifndef particle_h
#define particle_h

#include <iostream>
#include <memory>
#include <valarray>

#include "./scatterer.h"

namespace mc
{

class free_flight;

class particle
{

private:

  // pointer to scatter object for scattering the particle
  const scatterer* _scat_ptr=nullptr;

  // position of the particle
  std::valarray<double> _pos;

  // position of the particle in the previous time step, this is used for boundary collision detection
  std::valarray<double> _old_pos;

  // free flight time until next scattering event
  double _ff_time;

  // boolean determining the direction of the movement of particle in the cnt
  bool _heading_right;
  
  scatterer* _next_scat=nullptr;

  // 
  double _velocity;

public:
  particle() : _scat_ptr(nullptr), _pos({0, 0, 0}), _old_pos({0, 0, 0}), _ff_time(0), _heading_right(true), _velocity(0){};

  particle(const std::valarray<double>& pos, const scatterer* s, const double& velocity)
      : _scat_ptr(s), _pos(pos), _old_pos(pos), _velocity(velocity) {
    _ff_time = scat_ptr()->ff_time();
    _heading_right = std::rand()%2;
  };

  // perform free flight within the simulation domain
  void fly(double dt, const std::vector<scatterer>& s_list);

  // set the pointer to the scatterer object
  void set_scatterer(const scatterer* s) { _scat_ptr = s; };

  // return the pointer to the scatterer object
  const scatterer* scat_ptr() const { return _scat_ptr; };

  // get position of the particle
  const std::valarray<double>& pos() const { return _pos; };

  // get position of the particle
  const double& pos(const double& i) const { return _pos[i]; };

  // get old position of the particle
  const std::valarray<double>& old_pos() const { return _old_pos; };

  // get old position of the particle
  const double& old_pos(const int& i) const { return _old_pos[i]; };

  // set position of the particle and set the old position into _old_pos
  void set_pos(const std::valarray<double>& pos) { _pos = pos; };

  // // set an element of particle position and set the old position into _old_pos
  void set_pos(const int& i, const double& value) { _pos[i] = value; };

  // // set the old position of the particle.
  void set_old_pos(const std::valarray<double>& old_pos) { _old_pos = old_pos; };

  // return the free flight time until the next scattering
  const double& ff_time() const { return _ff_time; };

  // return the free flight time until the next scattering
  void set_ff_time(const double& value) { _ff_time = value; };

  // update the _ff_time by calling the underlying scatterer
  void get_ff_time() { _ff_time = _scat_ptr->ff_time(); };

  // step particle state for dt in time
  void step(double dt, const std::vector<scatterer>& s_list, const double& max_hop_radius);

}; //particle class

} //mc namespace

#endif // particle_h
