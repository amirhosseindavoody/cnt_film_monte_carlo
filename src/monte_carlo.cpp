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


	mc::t_float eff_mass = mc::elec_mass;
	std::shared_ptr<mc::free_flight> pilot = std::make_shared<mc::free_flight>();
	
	for (int i=0; i<num_particles; ++i)
	{
		// get random energy with correct distribution
		mc::t_float energy = -(3./2./_beta)*std::log(mc::get_rand_include_zero<mc::t_float>());
		mc::t_float mag = std::sqrt(energy*2.*eff_mass);
		

		// get uniformly distribution direction
		mc::t_float theta = std::acos(1.-2.*mc::get_rand_include_zero<mc::t_float>());
		mc::t_float phi = 2.*mc::pi*mc::get_rand_include_zero<mc::t_float>();

		mc::arr1d pos = {0., 0., 0.};
		mc::arr1d momentum = {mag*std::sin(theta)*std::cos(phi), mag*std::sin(theta)*std::sin(phi), mag*std::cos(theta)};
		_particles.emplace_back(pos, momentum, eff_mass, pilot);

	}
	_num_particles = _particles.size();
	_time = 0.;

};

// set the output directory and the output file name
void monte_carlo::process_command_line_args(int argc, char* argv[])
{
	namespace fs = std::experimental::filesystem;

	std::cout << "current path is " << fs::current_path() << std::endl;
	
	if (argc <= 1)
	{
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

	if (argc <= 2)
	{
		_output_file_path = _output_directory.path() / "output.dat";
	}
	else
	{
		_output_file_path = _output_directory.path() / std::string(argv[2]);
	}
	std::cout << "output file name:" << _output_file_path << std::endl;
};

// returns the number of particles
unsigned long int monte_carlo::num_particles()
{
	return _num_particles;
};

// step the simulation in time
void monte_carlo::step(mc::t_float dt)
{

	mc::t_float new_dt;
	for (std::list<mc::particle>::iterator it = _particles.begin(); it != _particles.end(); ++it)
	{
		new_dt = dt;
		
		while(it->get_ff_time() <= new_dt)
		{
			new_dt -= it->get_ff_time();
			it->fly(it->get_ff_time());
			it->scatter();
			it->update_ff_time();
		}

		it->fly(new_dt);
		it->get_ff_time() -= new_dt;
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
	for (std::list<mc::particle>::iterator it = _particles.begin(); it != _particles.end(); ++it)
	{
		file << (*it);
		file << "; ";
	}
	file << std::endl;
};

// get the mc simulation time
mc::t_float& monte_carlo::time()
{
	return _time;
};


} // namespace mc