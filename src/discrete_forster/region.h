#ifndef region_h
#define region_h

#include <iostream>
#include <array>
#include <memory>
#include <list>

#include "../helper/utility.h"
#include "./particle.h"
#include "./scatterer.h"

namespace mc
{

class region_class
{
private:
  unsigned _id; // this is a unique id for each region which is used for hashing
  arma::vec _lower_corner; // coordinate of the lower corner of the region
  arma::vec _upper_corner; // coordinate of the upper corner or the region
  double _volume; // volume of the region
  unsigned _number_of_scatterers;

  int _particle_flow_log; // this is the net number of particles flowing in (positive) or out (negative) of the region, the first component is the particle flow, the second number is the history.

  std::list<std::unique_ptr<particle>> _particles; // list of particles in the region
  std::list<std::unique_ptr<particle>> _new_particles; // list of particles newly entered the region
  std::vector<scatterer*> _scatterer_vector; // list of scatterers

  std::list<free_flight> _pilot_list; // list of free_flight objects

public:
  // default constructor
  region_class() {};

  // set boarders of the region
  void set_borders(const arma::vec& lower_corner, const arma::vec& upper_corner) {
    if (arma::any(lower_corner > upper_corner)) {
      throw std::out_of_range("invalid region definition: lower corner is larger than the upper corner!");
    }

  	_lower_corner = lower_corner;
  	_upper_corner = upper_corner;

  	_volume = arma::prod(_upper_corner-_lower_corner);

  };

  // decrease the net particle flow by one
  void loose_particle() {
    _particle_flow_log -= 1;
  };

  // increase the net particle flow by one
  void get_particle() {
    _particle_flow_log += 1;
  };

  // read the net particle flow
  const int& particle_flow() const {
    return _particle_flow_log;
  };

  // reset the particle flow counter to zero.
  void reset_particle_flow() {
    _particle_flow_log = 0;
  };

  // add a particle to the _particles list if the particle is in the region region
  typedef std::list<std::unique_ptr<particle>>::iterator pIterator;
  bool enlist(pIterator& particle_iterator, region_class* other_region) {
  	bool is_in_region = in_region(*(*particle_iterator));
  	if (is_in_region) {
  		auto prev_iterator = std::prev(particle_iterator,1); // get the iterator of the previous particle from the other region particle list
  		_new_particles.splice(_new_particles.end(), other_region->particles(), particle_iterator);
  		particle_iterator = prev_iterator; // now the particle iterator is the previous particle from the other region particle list
      other_region->loose_particle();
      get_particle();
  	}
  	return is_in_region;
  };

  // checks if a coordinate is inside the region
  bool in_region(const arma::vec& pos) {
    return (arma::all(pos>=_lower_corner) && arma::all(pos<=_upper_corner));
  };

  // checks if a particle is inside the region
  bool in_region(const mc::particle& p) {
    return (arma::all(p.pos()>=_lower_corner) && arma::all(p.pos()<=_upper_corner));
  };

  // gives the number of particles
  unsigned number_of_particles() const {
    return _particles.size() + _new_particles.size();
  };

  // gives an element of the _lower_corner
  const double& lower_corner(unsigned i) const {
    return _lower_corner(i);
  };

  // gives _lower_corner
  const arma::vec& lower_corner() const {
    return _lower_corner;
  };


  // gives an element of the _upper_corner
  const double& upper_corner(unsigned i) const {
    return _upper_corner[i];
  };

  // gives _upper_corner
  const arma::vec& upper_corner() const {
    return _upper_corner;
  };

  // populate the region with a certain number of particles
  void populate(const unsigned& number_of_particles) {
    // create the new particles by updating the previous particles and adding new ones or deleting the excess ones
    dump_new_particles();

    unsigned dice;

  	unsigned count=0;
  	auto p = _particles.begin();
  	while(count < number_of_particles) {
      dice = std::rand()%_number_of_scatterers;
  		if (p!= _particles.end()) {
        (*p)->set_pos(_scatterer_vector[dice]->pos());
        (*p)->set_old_pos(_scatterer_vector[dice]->pos());
        (*p)->set_scatterer(_scatterer_vector[dice]);
        (*p)->get_ff_time();
  			++p;
  		} else {
        _particles.push_back(std::make_unique<particle>( _scatterer_vector[dice]->pos(), &(_pilot_list.back()), _scatterer_vector[dice]));
  		}
  		count++;
  	}

  	_particles.erase(p, _particles.end());

  };

  // dump _new_particles into _particles list
  void dump_new_particles() {
  	_particles.splice(_particles.end(),_new_particles);
  };


  // return _particles list
  typedef std::list<std::unique_ptr<particle>> plist;
  plist& particles() {
    return _particles;
  };

  // get volume of the region
  const double& volume() const {
    return _volume;
  };

  // create a list of all scatterers that are inside this region
  typedef std::vector<scatterer> slist;
  void create_scatterer_vector(slist& all_scat_list) {
    for (auto& scat : all_scat_list) {
      if (in_region(scat.pos())) {
        _scatterer_vector.push_back(&scat);
      }
    }
    _number_of_scatterers = _scatterer_vector.size();
  };

  // get number of scatterers in the region
  unsigned number_of_scatterers() {
    return unsigned(_number_of_scatterers);
  }

  // create a list of all free_flight objects
  void create_pilot_list() {
    _pilot_list.push_back(free_flight());
  };

}; //region_class

} //mc namespace

#endif // region_h
