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

	_pilot = std::make_shared<mc::free_flight>();
	_scatterer = std::make_shared<mc::scatter>();

	_contacts[0].populate(_beta, 1100, _pilot, _scatterer);
	_contacts[1].populate(_beta, 100, _pilot, _scatterer);
	_bulk.populate(_beta,0, _pilot, _scatterer);

	_number_of_profile_sections = 10;
	_history_of_population_profiler = 0;
	_population_profile = std::vector<mc::t_uint>(_number_of_profile_sections, 0);
	_current_profile = std::vector<mc::t_float>(_number_of_profile_sections, 0.);
	_current[_contacts[0]] = 0;
	_current[_contacts[1]] = 0;
	_current[_bulk] = 0;

};

// set the output directory and the output file name
void monte_carlo::process_command_line_args(int argc, char* argv[])
{
	namespace fs = std::experimental::filesystem;

	std::cout << "current path is " << fs::current_path() << std::endl;

	if (argc <= 1)
	{
		// _output_directory.assign("/Users/amirhossein/research/test");
		_output_directory.assign("/home/amirhossein/research/test");
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

// step the simulation in time
void monte_carlo::step(mc::t_float dt)
{
	for (auto& contact: _contacts)
	{
		for (auto&& it=contact.particles().begin(); it!=contact.particles().end(); ++it)
		{
			it->step(dt,_domain);
			bool enlisted = _bulk.enlist(it,contact.particles());
			if (enlisted)
			{
				_current[contact] -= 1;
				_current[_bulk] += 1;
			}
		}
	}

	// move bulk to left_contact and right_contact
	for (auto&& it=_bulk.particles().begin(); it!=_bulk.particles().end(); ++it)
	{
		it->step(dt,_domain);
		for (auto& contact: _contacts)
		{
			bool enlisted = contact.enlist(it,_bulk.particles());
			if (enlisted)
			{
				_current[_bulk] -= 1;
				_current[contact] += 1;
				break;
			}
		}
	}

	// dump the new particles into the particles list in each region
	_bulk.dump_new_particles();
	for (auto& contact: _contacts)
	{
		contact.dump_new_particles();
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

	file << _time << " , " << number_of_particles() << " ; ";
	for (const auto& m_particle: _bulk.particles())
	{
		file << m_particle;
		file << "; ";
	}
	file << std::endl;
};

// calculate and save the population profile
void monte_carlo::update_profile(const mc::t_uint& max_history, std::fstream &population_file, std::fstream& current_file)
{
	if (_history_of_population_profiler < max_history)
	{
		_history_of_population_profiler += 1;

		mc::t_float length = (_bulk.upper_corner(2)-_bulk.lower_corner(2))/mc::t_float(_population_profile.size());
		for (const auto& m_particle: _bulk.particles())
		{
			int i = std::floor((m_particle.pos(2)-_bulk.lower_corner(2))/length);
			_population_profile[i] += 1;
			_current_profile[i] += m_particle.velocity(2);
		}
	}
	else
	{
		mc::t_float section_volume = _bulk.volume()/mc::t_float(_population_profile.size());

		if (! population_file.is_open())
		{
			population_file.open(_output_directory.path() / "population_profile.dat", std::ios::out);
			population_file << std::showpos << std::scientific;
		}
		population_file << time() << " ";
		for (auto& element: _population_profile)
		{
			population_file << mc::t_float(element)/mc::t_float(_history_of_population_profiler)/section_volume;
			population_file << " ";
			element = 0;
		}
		population_file << std::endl;

		if (! current_file.is_open())
		{
			current_file.open(_output_directory.path() / "current_profile.dat", std::ios::out);
			current_file << std::showpos << std::scientific;
		}
		current_file << time() << " ";
		for (auto& element: _current_profile)
		{
			current_file << mc::q0*element/mc::t_float(_history_of_population_profiler)/section_volume;
			current_file << " ";
			element = 0;
		}
		current_file << std::endl;

		_history_of_population_profiler = 0;
	}
};

} // namespace mc
