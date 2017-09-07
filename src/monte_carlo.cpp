#include <iostream>
#include <memory>
#include <experimental/filesystem>
#include <ctime>
#include <cmath>
#include <fstream>

#include "monte_carlo.h"
#include "particle.h"
#include "utility.h"

namespace mc
{

// constructor
monte_carlo::monte_carlo(unsigned long int num_particles):
	_left_contact({0,0,0}, {100.e-9, 100.e-9, 100.e-9}, 100),
	_bulk({0,0,100.e-9}, {100.e-9,100.e-9,900.e-9}, 0),
	_right_contact({0,0,900.e-9}, {100.e-9,100.e-9,1000.e-9})
{
	mc::init_random_number_generator();

	_temperature = 300;
	_beta = 1./(mc::kB*_temperature);
	_volume = {100.e-9, 100.e-9, 1000.e-9};
	_time = 0.;


	mc::t_float eff_mass = mc::elec_mass;
	mc::arr1d acceleration = {0., 0., 0.};
	std::shared_ptr<mc::free_flight> pilot = std::make_shared<mc::free_flight>(acceleration);
	std::shared_ptr<mc::scatter> scatterer = std::make_shared<mc::scatter>();

	mc::t_int id = 0;

	_particles.reserve(num_particles);

	for (int i=0; i<num_particles; ++i)
	{

		mc::arr1d pos;
		for (int j=0; j<pos.size(); ++j)
			pos[j] = _volume[j]*mc::get_rand_include_zero<mc::t_float>();

		// get random energy with correct distribution
		mc::t_float energy = -(3./2./_beta)*std::log(mc::get_rand_include_zero<mc::t_float>());
		mc::t_float velocity_magnitude = std::sqrt(energy*2./eff_mass);
		// get uniformly distribution direction
		mc::t_float theta = std::acos(1.-2.*mc::get_rand_include_zero<mc::t_float>());
		mc::t_float phi = 2.*mc::pi*mc::get_rand_include_zero<mc::t_float>();
		mc::arr1d velocity = {velocity_magnitude*std::sin(theta)*std::cos(phi), velocity_magnitude*std::sin(theta)*std::sin(phi), velocity_magnitude*std::cos(theta)};

		_particles.emplace_back(pos, velocity, eff_mass, pilot, scatterer, id);
		_bulk_particles_index.emplace_back(i);
		id ++;

	}
	_num_particles = _particles.size();

};

// set the output directory and the output file name
void monte_carlo::process_command_line_args(int argc, char* argv[])
{
	namespace fs = std::experimental::filesystem;

	std::cout << "current path is " << fs::current_path() << std::endl;

	if (argc <= 1)
	{
		_output_directory.assign("/Users/amirhossein/research/test");
		// _output_directory.assign("/home/amirhossein/research/test");
	}
	else
	{
		_output_directory.assign(argv[1]);
	}

	if (not fs::exists(_output_directory.path()))
	{
		std::cout << "warning: output directory does NOT exist!!!" << std::endl;
		std::cout << "output directory: " << _output_directory.path() << std::endl;
		fs::create_directories(_output_directory.path());
		// std::exit(EXIT_FAILURE);
	}

	if (not fs::is_directory(_output_directory.path()))
	{
		std::cout << "error: output path is NOT a directory!!!" << std::endl;
		std::cout << "output path: " << _output_directory.path() << std::endl;
		std::exit(EXIT_FAILURE);
	}

	if (not fs::is_empty(_output_directory.path()))
	{
		std::cout << "warning: output directory is NOT empty!!!" << std::endl;
		std::cout << "output directory: " << _output_directory.path() << std::endl;
		std::cout << "deleting the existing directory!!!" << std::endl;
		fs::remove_all(_output_directory.path());
		fs::create_directories(_output_directory.path());
		// std::exit(EXIT_FAILURE);
	}

};

// returns the number of particles
unsigned long int monte_carlo::num_particles()
{
	return _num_particles;
};

// step the simulation in time
void monte_carlo::step(mc::t_float dt)
{
	// step particles in the bulk section
	for (const auto& index: _bulk_particles_index)
	{
		_particles[index].step(dt,_volume);
	}

	// step particles in the right contact section
	for (const auto& index: _right_contact_particles_index)
	{
		_particles[index].step(dt,_volume);
	}

	//step particles in the left contact sections
	for (const auto& index: _left_contact_particles_index)
	{
		_particles[index].step(dt,_volume);
	}

	_time += dt;

};

// write the state of the simulation in the output file
void monte_carlo::write_state(std::fstream &file)
{
	if (! file.is_open())
	{
		file.open(_output_directory.path() / "output.dat", std::ios::out);
		file << std::showpos << std::scientific;
	}

	file << _time << " , " << num_particles() << " ; ";
	for (const auto& m_particle: _particles)
	{
		file << m_particle;
		file << "; ";
	}
	file << std::endl;
};

// get the mc simulation time
mc::t_float& monte_carlo::time()
{
	return _time;
};

// update the list of particle indices for active, inactive, bulk, and contact regions.
void monte_carlo::update_particle_list()
{
	// move left_contact to bulk and right_contact
	for (auto it=_left_contact_particles_index.begin(); it!=_left_contact_particles_index.end(); ++it)
	{
		if (_bulk.in_region(_particles[*it].pos()))
		{
			--it;
			_bulk_particles_index.splice(_bulk_particles_index.end(),_left_contact_particles_index, std::next(it,1));
		}
		else if (_right_contact.in_region(_particles[*it].pos()))
		{
			--it;
			_right_contact_particles_index.splice(_right_contact_particles_index.end(),_left_contact_particles_index, std::next(it,1));
		}
	}

	// move bulk to left_contact and right_contact
	for (auto it=_bulk_particles_index.begin(); it!=_bulk_particles_index.end(); ++it)
	{
		if (_left_contact.in_region(_particles[*it].pos()))
		{
			--it;
			_left_contact_particles_index.splice(_left_contact_particles_index.end(),_bulk_particles_index, std::next(it,1));
		}
		else if (_right_contact.in_region(_particles[*it].pos()))
		{
			--it;
			_right_contact_particles_index.splice(_right_contact_particles_index.end(),_bulk_particles_index, std::next(it,1));
		}
	}

	// move right_contact to bulk and left_contact
	for (auto it=_right_contact_particles_index.begin(); it!=_right_contact_particles_index.end(); ++it)
	{
		if (_left_contact.in_region(_particles[*it].pos()))
		{
			--it;
			_left_contact_particles_index.splice(_left_contact_particles_index.end(),_right_contact_particles_index, std::next(it,1));
		}
		else if (_bulk.in_region(_particles[*it].pos()))
		{
			--it;
			_bulk_particles_index.splice(_bulk_particles_index.end(),_right_contact_particles_index, std::next(it,1));
		}
	}

	// // sanity check to see if we did everything correctly
	// bool correct_particle_region = true;
	// for (auto it=_left_contact_particles_index.begin(); it!=_left_contact_particles_index.end(); ++it)
	// {
	// 	if (! _left_contact.in_region(_particles[*it].pos()))
	// 	{
	// 		std::cout << "particle is supposed to be in the left contact: id= " << _particles[*it].id() << std::endl;
	// 		correct_particle_region = false;
	// 	}
	// }
	// for (auto it=_right_contact_particles_index.begin(); it!=_right_contact_particles_index.end(); ++it)
	// {
	// 	if (! _right_contact.in_region(_particles[*it].pos()))
	// 	{
	// 		std::cout << "particle is supposed to be in the right contact: id= " << _particles[*it].id() << std::endl;
	// 		correct_particle_region = false;
	// 	}
	// }
	// for (auto it=_bulk_particles_index.begin(); it!=_bulk_particles_index.end(); ++it)
	// {
	// 	if (! _bulk.in_region(_particles[*it].pos()))
	// 	{
	// 		std::cout << "particle is supposed to be in the bulk: id= " << _particles[*it].id() << std::endl;
	// 		correct_particle_region = false;
	// 	}
	// }
	//
	// if (! correct_particle_region)
	// {
	// 	std::cin.ignore();
	// }

};

} // namespace mc
