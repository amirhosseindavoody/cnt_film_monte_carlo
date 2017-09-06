#include <iostream>

#include "contact.h"

namespace mc
{

// constructor
contact::contact(const mc::arr1d lower_corner, const mc::arr1d upper_corner, const mc::t_int number_of_expected_particles)
{
	_lower_corner = lower_corner;
	_upper_corner = upper_corner;

	for (int i=0; i<_lower_corner.size(); ++i)
	{
		if (_lower_corner[i] >= _upper_corner[i])
		{
			std::cout << "ERROR: invalid contact definition--> lower corner is larger than the upper corner!" << std::endl;
			exit(0);
		}
	}

	_volume = 1.0;
	for (int i=0; i<_lower_corner.size(); ++i)
	{
		_volume = _volume*(_upper_corner[i]-_lower_corner[i]);
	}

	_number_of_expected_particles = number_of_expected_particles;
	_number_of_actual_particles = 0;
};

// add a particle to the _particles list if the particle is in the contact region
void contact::enlist(const std::list<mc::particle>::iterator particle_ptr)
{
	bool is_in_contact = true;

	for (int i=0; i<_lower_corner.size(); ++i)
	{
		if ((particle_ptr->pos(i) < _lower_corner[i]) || (particle_ptr->pos(i) > _upper_corner[i]))
		{
			is_in_contact = false;
		}
	}

	if (is_in_contact)
	{
		std::cout << "adding particle to contact area list!!!" << std::endl;
		_particle_iterators.emplace_back(particle_ptr);
		std::cin.ignore();
	}

};

} // mc namespace
