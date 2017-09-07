#include <iostream>

#include "region.h"

namespace mc
{

// constructor
region::region(const mc::arr1d lower_corner, const mc::arr1d upper_corner, const mc::t_int number_of_expected_particles)
{
	_lower_corner = lower_corner;
	_upper_corner = upper_corner;

	for (int i=0; i<_lower_corner.size(); ++i)
	{
		if (_lower_corner[i] >= _upper_corner[i])
		{
			std::cout << "ERROR: invalid region definition--> lower corner is larger than the upper corner!" << std::endl;
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

// add a particle to the _particles list if the particle is in the region region
void region::enlist(const std::list<mc::particle>::iterator particle_ptr)
{
	bool is_in_region = true;

	for (int i=0; i<_lower_corner.size(); ++i)
	{
		if ((particle_ptr->pos(i) < _lower_corner[i]) || (particle_ptr->pos(i) > _upper_corner[i]))
		{
			is_in_region = false;
		}
	}

	if (is_in_region)
	{
		std::cout << "adding particle to region area list!!!" << std::endl;
		_particle_iterators.emplace_back(particle_ptr);
		std::cin.ignore();
	}

};

// checks if a coordinate is inside the region
bool region::in_region(const mc::arr1d& pos)
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

} // mc namespace
