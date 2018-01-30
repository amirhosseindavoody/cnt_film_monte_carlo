#ifndef discrete_forster_region_h
#define discrete_forster_region_h

#include <iostream>
#include <array>
#include <memory>
#include <list>

#include "../helper/utility.h"
#include "discrete_forster_particle.h"

namespace mc
{

class discrete_forster_region
{
public:
  typedef mc::discrete_forster_particle t_particle; // particle type
  typedef mc::discrete_forster_free_flight t_ff; // free_flight type
  typedef mc::discrete_forster_region t_region; // region type
  typedef mc::discrete_forster_scatter t_scatter; // scatter

private:
  mc::t_uint _id; // this is a unique id for each region which is used for hashing
  mc::arr1d _lower_corner; // coordinate of the lower corner of the region
  mc::arr1d _upper_corner; // coordinate of the upper corner or the region
  mc::t_float _volume; // volume of the region
  mc::t_float _number_of_scatterers;

  mc::t_int _particle_flow_log; // this is the net number of particles flowing in (positive) or out (negative) of the region, the first component is the particle flow, the second number is the history.

  std::list<std::unique_ptr<t_particle>> _particles; // list of particles in the region
  std::list<std::unique_ptr<t_particle>> _new_particles; // list of particles newly entered the region
  std::vector<std::shared_ptr<t_scatter>> _scatterer_vector; // list of scatterers
  std::list<std::shared_ptr<t_ff>> _pilot_list; // list of free_flight objects

public:
  void set_borders(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner) // set boarders of the region
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
  void loose_particle() // decrease the net particle flow by one
  {
    _particle_flow_log -= 1;
  };
  void get_particle() // increase the net particle flow by one
  {
    _particle_flow_log += 1;
  };
  const mc::t_int& particle_flow() const // read the net particle flow
  {
    return _particle_flow_log;
  };
  void reset_particle_flow() // reset the particle flow counter to zero.
  {
    _particle_flow_log = 0;
  }
  bool enlist(std::list<std::unique_ptr<t_particle>>::iterator& particle_iterator, t_region& other_region) // add a particle to the _particles list if the particle is in the region region
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
  bool in_region(const mc::arr1d& pos) // checks if a coordinate is inside the region
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
  bool in_region(const t_particle& p) // checks if a particle is inside the region
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
  mc::t_uint number_of_particles() const // gives the number of particles
  {
    return _particles.size() + _new_particles.size();
  };
  const mc::t_float& lower_corner(mc::t_uint i) const // gives an element of the _lower_corner
  {
    return _lower_corner[i];
  };
  const mc::arr1d& lower_corner() const // gives _lower_corner
  {
    return _lower_corner;
  };
  const mc::t_float& upper_corner(mc::t_uint i) const // gives an element of the _upper_corner
  {
    return _upper_corner[i];
  };
  const mc::arr1d& upper_corner() const // gives _upper_corner
  {
    return _upper_corner;
  };
  void populate(const mc::t_uint& number_of_particles) // populate the region with a certain number of particles
  {
    // create the new particles by updating the previous particles and adding new ones or deleting the excess ones
    dump_new_particles();

    mc::t_uint dice;

  	mc::t_uint count=0;
  	auto it = _particles.begin();
  	while(count < number_of_particles)
  	{
      dice = (_number_of_scatterers*mc::get_rand_include_zero<mc::t_float>());
  		if (it!= _particles.end())
  		{
        (*it)->set_pos(_scatterer_vector[dice]->pos());
        (*it)->set_scatterer(_scatterer_vector[dice]);
        (*it)->get_ff_time();
  			it++;
  		}
  		else
  		{
        _particles.push_back(std::make_unique<t_particle>( _scatterer_vector[dice]->pos(), _pilot_list.back(), _scatterer_vector[dice]));
  		}
  		count++;
  	}

  	_particles.erase(it, _particles.end());

  };
  void dump_new_particles() // dump _new_particles into _particles list
  {
  	_particles.splice(_particles.end(),_new_particles);
  };
  std::list<std::unique_ptr<t_particle>>& particles() // return _particles list
  {
  	return _particles;
  };
  const mc::t_float& volume() const // get volume of the region
  {
    return _volume;
  };
  void create_scatterer_vector(const std::list<std::shared_ptr<t_scatter>>& all_scat_list) // create a list of all scatterers that are inside this region
  {
    for (const auto& scat : all_scat_list)
    {
      if (in_region(scat->pos()))
      {
        _scatterer_vector.push_back(scat);
      }
    }
    _number_of_scatterers = _scatterer_vector.size();
    // std::cout << "scatter list size: " << _scatterer_vector.size() << std::endl;
  }
  void create_pilot_list() // create a list of all free_flight objects
  {
    _pilot_list.push_back(std::make_shared<t_ff>());
  }

}; //discrete_forster_region class

} //mc namespace

#endif // discrete_forster_region_h
