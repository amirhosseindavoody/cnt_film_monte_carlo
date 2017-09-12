#include <iostream>

#include "region.h"

namespace mc
{

// constructor
region::region(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner)
{
	set_borders(lower_corner, upper_corner);
};

// add a particle to the _particles list if the particle is in the region region
bool region::enlist(std::list<mc::particle>::iterator& particle_iterator, std::list<mc::particle>& other_region_particles)
{
	bool is_in_region = in_region(*particle_iterator);
	if (is_in_region)
	{
		auto prev_iterator = std::prev(particle_iterator,1); // get the iterator of the previous particle from the other region particle list
		_new_particles.splice(_new_particles.end(), other_region_particles, particle_iterator);
		particle_iterator = prev_iterator; // now the particle iterator is the previous particle from the other region particle list
	}

	return is_in_region;
};

// populate the region with a certain number of particles
void region::populate(const mc::t_float& beta, const mc::t_uint& number_of_particles, const std::shared_ptr<mc::free_flight>& pilot, const std::shared_ptr<mc::scatter>& scatterer)
{
	_number_of_expected_particles = number_of_particles;

	mc::t_float eff_mass = mc::elec_mass;
	mc::arr1d acceleration = {0., 0., 0.};
	// std::shared_ptr<mc::free_flight> pilot = std::make_shared<mc::free_flight>(acceleration);
	// std::shared_ptr<mc::scatter> scatterer = std::make_shared<mc::scatter>();

	mc::t_int id = 0;

	// _particles.clear();
	// _new_particles.clear();

	for (int i=0; i<_number_of_expected_particles; ++i)
	{
		mc::arr1d pos;
		for (int j=0; j<pos.size(); ++j)
			pos[j] = _lower_corner[j]+(_upper_corner[j]-_lower_corner[j])*mc::get_rand_include_zero<mc::t_float>();

		// get random energy with correct distribution
		mc::t_float energy = -(3./2./beta)*std::log(mc::get_rand_include_zero<mc::t_float>());
		mc::t_float velocity_magnitude = std::sqrt(energy*2./eff_mass);
		// get uniformly distribution direction
		mc::t_float theta = std::acos(1.-2.*mc::get_rand_include_zero<mc::t_float>());
		mc::t_float phi = 2.*mc::pi*mc::get_rand_include_zero<mc::t_float>();
		mc::arr1d velocity = {velocity_magnitude*std::sin(theta)*std::cos(phi), velocity_magnitude*std::sin(theta)*std::sin(phi), velocity_magnitude*std::cos(theta)};

		_particles.emplace_back(pos, velocity, eff_mass, pilot, scatterer, id);
		++id;
	}

	dump_new_particles();

	mc::t_uint count=0;
	auto it = _particles.begin();
	while((count < _number_of_expected_particles) and (it != _particles.end()))
	{
		mc::arr1d pos;
		for (int j=0; j<pos.size(); ++j)
			pos[j] = _lower_corner[j]+(_upper_corner[j]-_lower_corner[j])*mc::get_rand_include_zero<mc::t_float>();

		// get random energy with correct distribution
		mc::t_float energy = -(3./2./beta)*std::log(mc::get_rand_include_zero<mc::t_float>());
		mc::t_float velocity_magnitude = std::sqrt(energy*2./eff_mass);
		// get uniformly distribution direction
		mc::t_float theta = std::acos(1.-2.*mc::get_rand_include_zero<mc::t_float>());
		mc::t_float phi = 2.*mc::pi*mc::get_rand_include_zero<mc::t_float>();
		mc::arr1d velocity = {velocity_magnitude*std::sin(theta)*std::cos(phi), velocity_magnitude*std::sin(theta)*std::sin(phi), velocity_magnitude*std::cos(theta)};

		// _particles.emplace_back(pos, velocity, eff_mass, pilot, scatterer, id);
		it->set_pos(pos);
		it->set_velocity(velocity);
		it->kin_energy();

		++id;

		it ++;
		count ++;
	}

};

// dump _new_particles into _particles list
void region::dump_new_particles()
{
	_particles.splice(_particles.end(),_new_particles);
};

// return _particles list
std::list<mc::particle>& region::particles()
{
	return _particles;
};

// set boarders of the region
void region::set_borders(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner)
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

} // mc namespace
