#include <iostream>

#include "contact.h"

namespace mc
{

// constructor
contact::contact(const mc::arr1d lower_corner, const mc::arr1d upper_corner)
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
};

} // mc namespace
