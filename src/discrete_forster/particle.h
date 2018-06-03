#ifndef particle_h
#define particle_h

#include <iostream>
#include <array>
#include <memory>
#include <armadillo>

#include "./free_flight.h"
#include "./scatterer.h"

namespace mc
{

class free_flight;

class particle
{

private:

  // pointer to free_flight object for driving the particle
  const free_flight* _pilot=nullptr;

  // pointer to scatter object for scattering the particle
  const scatterer* _scat_ptr=nullptr;

  // position of the particle
  arma::vec _pos;

  // position of the particle in the previous time step, this is used for boundary collision detection
  arma::vec _old_pos;

  // free flight time until next scattering event
  double _ff_time;

public:
  particle() : _pilot(nullptr), _scat_ptr(nullptr), _pos({0, 0, 0}), _old_pos({0, 0, 0}), _ff_time(0){};

  particle(const arma::vec& pos, const free_flight* pilot, const scatterer* s)
      : _pilot(pilot), _scat_ptr(s), _pos(pos), _old_pos(pos) {
    _ff_time = scat_ptr()->ff_time();
  };

  // perform free flight within the simulation domain
  void fly(const double& dt, const std::pair<arma::vec, arma::vec>& domain);

  // set the pilot free_flight pointer object
  void set_pilot(free_flight* pilot) { _pilot = pilot; };

  // get the pilot free_flight pointer
  const free_flight* pilot() const { return _pilot; };

  // set the pointer to the scatterer object
  void set_scatterer(const scatterer* s) { _scat_ptr = s; };

  // return the pointer to the scatterer object
  const scatterer* scat_ptr() const { return _scat_ptr; };

  // get position of the particle
  const arma::vec& pos() const { return _pos; };

  // get position of the particle
  const double& pos(const double& i) const { return _pos(i); };

  // get old position of the particle
  const arma::vec& old_pos() const { return _old_pos; };

  // get old position of the particle
  const double& old_pos(const int& i) const { return _old_pos(i); };

  // set position of the particle and set the old position into _old_pos
  void set_pos(const arma::vec& pos) { _pos = pos; };

  // // set an element of particle position and set the old position into _old_pos
  void set_pos(const int& i, const double& value) { _pos(i) = value; };

  // // set the old position of the particle.
  void set_old_pos(const arma::vec& old_pos) { _old_pos = old_pos; };

  // return the free flight time until the next scattering
  const double& ff_time() const { return _ff_time; };

  // return the free flight time until the next scattering
  void set_ff_time(const double& value) { _ff_time = value; };

  // update the _ff_time by calling the underlying scatterer
  void get_ff_time() { _ff_time = _scat_ptr->ff_time(); };

  // step particle state for dt in time
  void step(double dt,
            const std::pair<arma::vec, arma::vec>& domain,
            const double& max_hop_radius);

}; //particle class

} //mc namespace

#endif // particle_h
