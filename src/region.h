#ifndef region_h
#define region_h

#include <iostream>
#include <array>
#include <memory>

#include "ff.h"
#include "utility.h"
#include "particle.h"
#include "gas_particle.h"
#include "forster_particle.h"

namespace mc
{

class region
{
private:
  mc::t_uint _id; // this is a unique id for each region which is used for hashing
  mc::arr1d _lower_corner; // coordinate of the lower corner of the region
  mc::arr1d _upper_corner; // coordinate of the upper corner or the region
  mc::t_float _volume; // volume of the region region

  mc::t_int _particle_flow_log; // this is the net number of particles flowing in (positive) or out (negative) of the region, the first component is the particle flow, the second number is the history.

  typedef std::unique_ptr<mc::particle> particle_ptr; // particle pointers
  typedef std::list<particle_ptr> particle_ptr_list; // list of particle pointers

  particle_ptr_list _particles; // list of particles in the region
  particle_ptr_list _new_particles; // list of particles newly entered the region

public:
  // region(const mc::t_int& id = static_cast<mc::t_uint>(std::rand()), const mc::arr1d& lower_corner = {0,0,0}, const mc::arr1d& upper_corner = {0,0,0}) // constructor
  region(const mc::arr1d& lower_corner = {0,0,0}, const mc::arr1d& upper_corner = {0,0,0}, const mc::t_int& id = std::rand()): // constructor
    _id(id)
  {
  	set_borders(lower_corner, upper_corner);
  };
  inline void set_borders(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner) // set boarders of the region
  {
  	_lower_corner = lower_corner;
  	_upper_corner = upper_corner;

  	for (int i=0; i<_lower_corner.size(); ++i)
  	{
  		if (_lower_corner[i] > _upper_corner[i])
  		{
  			std::cout << "\n***\nERROR: invalid region definition--> lower corner is larger than the upper corner!" << std::endl;
  			std::cout << " ,lower corner[" << i <<"]: " << _lower_corner[i] << "\n ,upper corner[" << i << "]: " << _upper_corner[i] << "\n***\n";
  			// exit(0);
  		}
  	}

  	_volume = 1.0;
  	for (int i=0; i<_lower_corner.size(); ++i)
  	{
  		_volume = _volume*(_upper_corner[i]-_lower_corner[i]);
  	}
  };
  inline void loose_particle() // decrease the net particle flow by one
  {
    _particle_flow_log -= 1;
  };
  inline void get_particle() // increase the net particle flow by one
  {
    _particle_flow_log += 1;
  };
  inline const mc::t_int& particle_flow() const // read the net particle flow
  {
    return _particle_flow_log;
  };
  inline void reset_particle_flow() // reset the particle flow counter to zero.
  {
    _particle_flow_log = 0;
  }
  inline bool enlist(particle_ptr_list::iterator& particle_iterator, mc::region& other_region) // add a particle to the _particles list if the particle is in the region region
  {
  	bool is_in_region = in_region(*(*particle_iterator));
  	if (is_in_region)
  	{
  		auto prev_iterator = std::prev(particle_iterator,1); // get the iterator of the previous particle from the other region particle list
  		_new_particles.splice(_new_particles.end(), other_region.particles(), particle_iterator);
  		particle_iterator = prev_iterator; // now the particle iterator is the previous particle from the other region particle list
      other_region.loose_particle();
      get_particle();
  	}
  	return is_in_region;
  };
  inline bool in_region(const mc::arr1d& pos) // checks if a coordinate is inside the region
  {
  	for (int i=0; i<_lower_corner.size(); ++i)
  	{
  		if (pos[i] < _lower_corner[i])
  		{
  			return false;
  		}
  		if (pos[i] > _upper_corner[i])
  		{
  			return false;
  		}
  	}
  	return true;
  };
  inline bool in_region(const mc::particle& p) // checks if a particle is inside the region
  {
  	for (int i=0; i<_lower_corner.size(); ++i)
  	{
  		if (p.pos(i) < _lower_corner[i])
  		{
  			return false;
  		}
  		if (p.pos(i) > _upper_corner[i])
  		{
  			return false;
  		}
  	}
  	return true;
  };
  inline mc::t_uint number_of_particles() // gives the number of particles
  {
    return _particles.size() + _new_particles.size();
  };
  inline const mc::t_float& lower_corner(mc::t_uint i) const // gives an element of the _lower_corner
  {
    return _lower_corner[i];
  };
  inline const mc::arr1d& lower_corner() const // gives _lower_corner
  {
    return _lower_corner;
  };
  inline const mc::t_float& upper_corner(mc::t_uint i) const // gives an element of the _upper_corner
  {
    return _upper_corner[i];
  };
  inline const mc::arr1d& upper_corner() const // gives _upper_corner
  {
    return _upper_corner;
  };
  inline void populate(const mc::t_float& beta, const mc::t_uint& number_of_particles, const std::shared_ptr<mc::free_flight>& pilot, const std::shared_ptr<mc::scatter>& scatterer) // populate the region with a certain number of particles
  {
    // create the new particles by updating the previous particles and adding new ones or deleting the excess ones
    dump_new_particles();
  	mc::t_int id = 0;


  	mc::t_uint count=0;
  	auto it = _particles.begin();
  	while(count < number_of_particles)
  	{
      mc::arr1d pos;
  		for (int j=0; j<pos.size(); ++j)
      {
        pos[j] = _lower_corner[j]+(_upper_corner[j]-_lower_corner[j])*mc::get_rand_include_zero<mc::t_float>();
      }
  		mc::arr1d velocity = mc::get_rand_velocity(beta, mc::elec_mass);

  		if (it!= _particles.end())
  		{
  			(*it)->set_pos(pos);
  			(*it)->set_velocity(velocity);
  			(*it)->kin_energy();
  			it++;
  		}
  		else
  		{
        // _particles.push_back(std::make_unique<mc::gas_particle>(pos, velocity, eff_mass, pilot, scatterer, id));
        _particles.push_back(std::make_unique<mc::forster_particle>(pos, mc::elec_mass, pilot, scatterer, id));
  		}

  		id++;
  		count++;
  	}

  	_particles.erase(it, _particles.end());

  };
  inline void dump_new_particles() // dump _new_particles into _particles list
  {
  	_particles.splice(_particles.end(),_new_particles);
  };
  inline particle_ptr_list& particles() // return _particles list
  {
  	return _particles;
  };
  inline const mc::t_float& volume() const // get volume of the region
  {
    return _volume;
  };
  inline bool operator < (const mc::region& rhs) const // comparision operator for comparing regions for use in map data structures
  {
    return (_id < rhs._id);
  };

}; //region class

} //mc namespace

#endif // region_h
