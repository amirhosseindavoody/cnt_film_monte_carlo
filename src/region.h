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
  mc::t_int _number_of_actual_particles; // actual number of particles in the region region
  std::list<std::list<mc::particle>::iterator> _particle_iterators; // list of pointers to the particles in the region region
  // std::set

public:
  region(const mc::arr1d lower_corner, const mc::arr1d upper_corner, const mc::t_int number_of_expected_particles = 0);
  void enlist(const std::list<mc::particle>::iterator particle_ptr); // add a particle to the _particles list if the particle is in the region region
  bool in_region(const mc::arr1d& pos); // checks if a coordinate is inside the region

}; //region class

} //mc namespace

#endif // region_h
