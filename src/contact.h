#ifndef contact_h
#define contact_h

#include <iostream>
#include <array>
#include <memory>

#include "utility.h"

namespace mc
{

class contact
{
  mc::arr1d _lower_corner; // coordinate of the lower corner of the contact
  mc::arr1d _upper_corner; // coordinate of the upper corner or the contact
  mc::t_float _volume; // volume of the contact region
  mc::t_int _number_of_expected_particles; // number of expected particles in the contact region

private:
  contact(const mc::arr1d lower_corner, const mc::arr1d upper_corner, const mc::t_int number_of_expected_particles);

public:

}; //contact class

} //mc namespace

#endif // contact_h
