#ifndef region_h
#define region_h

#include <iostream>
#include <array>
#include <memory>

#include "utility.h"
#include "particle.h"

namespace mc
{

class region
{
private:
  mc::arr1d _lower_corner; // coordinate of the lower corner of the region
  mc::arr1d _upper_corner; // coordinate of the upper corner or the region
  mc::t_float _volume; // volume of the region region
  mc::t_int _number_of_expected_particles; // number of expected particles in the region region
  std::list<mc::particle> _particles; // list of particles in the region
  std::list<mc::particle> _new_particles; // list of particles newly entered the region

public:
  region(const mc::arr1d& lower_corner = {0,0,0}, const mc::arr1d& upper_corner = {0,0,0}); // constructor
  void set_borders(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner); // set boarders of the region
  bool enlist(std::list<mc::particle>::iterator& particle_iterator, std::list<mc::particle>& other_region_particles); // add a particle to the _particles list if the particle is in the region region
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
    return _particles.size();
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
  void populate(const mc::t_float& beta, const mc::t_uint& number_of_particles, const std::shared_ptr<mc::free_flight>& pilot = std::make_shared<mc::free_flight>(), const std::shared_ptr<mc::scatter>& scatterer = std::make_shared<mc::scatter>()); // populate the region with a certain number of particles
  void dump_new_particles(); // dump _new_particles into _particles list
  std::list<mc::particle>& particles(); // return _particles list


}; //region class

} //mc namespace

#endif // region_h
