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
monte_carlo::monte_carlo(unsigned long int num_particles)
{
	mc::init_random_number_generator();

	_temperature = 300;
	_beta = 1./(mc::kB*_temperature);
	_time = 0.;

	_contacts.emplace_back(mc::arr1d({0,0,0}), mc::arr1d({100.e-9, 100.e-9, 100.e-9}));
	_contacts.emplace_back(mc::arr1d({0,0,900.e-9}), mc::arr1d({100.e-9,100.e-9,1000.e-9}));
	_bulk.set_borders({0,0,100.e-9}, {100.e-9,100.e-9,900.e-9});

	_domain.first = {0., 0., 0.};
	_domain.second = {100.e-9, 100.e-9, 1000.e-9};

	_contacts[0].populate(_beta, 1100);
	_contacts[1].populate(_beta,100);
	_bulk.populate(_beta,0);

	_population_profile = std::vector<mc::t_uint>(10, 0);
	_history_of_population_profiler = 0;

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
mc::t_uint monte_carlo::num_particles()
{
	mc::t_uint number_of_particles = _bulk.number_of_particles();
	for (auto& contact: _contacts)
	{
		number_of_particles += contact.number_of_particles();
	}
	return number_of_particles;
};

// step the simulation in time
void monte_carlo::step(mc::t_float dt)
{
	for (auto& contact: _contacts)
	{
		for (auto&& it=contact.particles().begin(); it!=contact.particles().end(); ++it)
		{
			it->step(dt,_domain);
			_bulk.enlist(it,contact.particles());
		}
	}

	// move bulk to left_contact and right_contact
	for (auto&& it=_bulk.particles().begin(); it!=_bulk.particles().end(); ++it)
	{
		it->step(dt,_domain);
		for (auto& contact: _contacts)
		{
			bool enlisted = contact.enlist(it,_bulk.particles());
			if (enlisted)	break;
		}
	}

	// dump the new particles into the particles list in each region
	// for (auto& contact: _contacts)
	// {
	// 	contact.dump_new_particles();
	// }
	_bulk.dump_new_particles();

	_contacts[0].populate(_beta, 1100);
	_contacts[1].populate(_beta,100);

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
	for (const auto& m_particle: _bulk.particles())
	{
		file << m_particle;
		file << "; ";
	}
	file << std::endl;
};

// calculate and save the population profile
void monte_carlo::update_population_profile(const mc::t_uint& max_history, std::fstream &file)
{
	if (_history_of_population_profiler < max_history)
	{
		_history_of_population_profiler += 1;

		mc::t_float length = (_bulk.upper_corner(2)-_bulk.lower_corner(2))/mc::t_float(_population_profile.size());

		for (const auto& m_particle: _bulk.particles())
		{
			int i = std::floor((m_particle.pos(2)-_bulk.lower_corner(2))/length);
			_population_profile[i] += 1;
		}
	}
	else
	{
		if (! file.is_open())
		{
			file.open(_output_directory.path() / "population_profile.dat", std::ios::out);
			file << std::showpos << std::scientific;
		}

		file << time() << " ";
		for (auto& element: _population_profile)
		{
			file << mc::t_float(element)/mc::t_float(_history_of_population_profiler);
			file << " ";
			element = 0;
		}
		file << std::endl;
		_history_of_population_profiler = 0;
	}
};

} // namespace mc
