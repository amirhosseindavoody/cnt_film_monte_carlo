#ifndef monte_carlo_h
#define monte_carlo_h

#include <iostream>
#include <list>
#include <experimental/filesystem>
#include <fstream>
#include <map>
#include <chrono>

#include "particle.h"
#include "utility.h"
#include "region.h"

namespace mc
{

class monte_carlo
{
private:

	// simulation parameters
	mc::t_uint _num_particles; // number of particles
	mc::t_float _time; // total simulation time that has elapsed
	mc::t_float _temperature; // temperature of the simulation
	mc::t_float _beta; // 1/kB*T is the inverse of the thermal energy
	std::pair<mc::arr1d, mc::arr1d> _domain; // the physical extent of the simulation domain

	std::shared_ptr<mc::free_flight> _pilot; // pilot for the free flight step for particles in the simulation
	std::shared_ptr<mc::scatter> _scatterer; // scatterer for all particles in the simulation

	// defining various domains in the simulation
	mc::region _bulk;
	std::vector<mc::region> _contacts;

	// output files properties
	std::experimental::filesystem::directory_entry _output_directory; // this is the address of the output_directory
	mc::t_uint _number_of_profile_sections; // number of profile sections used for profiling
	std::pair<mc::t_uint, std::vector<mc::t_uint>> _population_profile; // this is the population profile of particles through the simulation domain along the z-axis
	// mc::t_uint _history_of_population_profiler; // this is the number of steps that the _population_profile have been recorded for considering the running average
	std::pair<mc::t_uint, std::vector<mc::t_float>> _current_profile; // this is the average current profile along the simulation domain
	mc::t_uint _history_of_region_currents; // this is the number of steps that the net _current in the regions have been recorded
	// std::map<mc::region, mc::t_int> _current; // this is the net current density into each reagion

public:

	inline monte_carlo() // constructor
	{
		mc::init_random_number_generator();

		_temperature = 300;
		_beta = 1./(mc::kB*_temperature);
		_time = 0.;

		// set simulation geometry parameters, the rest are set automatically
		mc::t_float contact_length = 100.e-9;
		mc::t_float bulk_length = 1000.e-9;
		mc::t_float cross_section_1 = 100.e-9;
		mc::t_float cross_section_2 = 100.e-9;
		mc::t_float total_length = 2.*contact_length+bulk_length;

		// set simulation domain for boundary reflection
		_domain.first = {0., 0., 0.};
		_domain.second = {cross_section_1, cross_section_2, total_length};

		// upper and lower corners of the contact and bulk
		mc::arr1d contact_1_lower_corner = {0., 0., 0.};
		mc::arr1d contact_1_upper_corner = {cross_section_1, cross_section_2, contact_length};
		mc::arr1d bulk_lower_corner = {0., 0., contact_length};
		mc::arr1d bulk_upper_corner = {cross_section_1, cross_section_2, contact_length+bulk_length};
		mc::arr1d contact_2_lower_corner = {0., 0., contact_length+bulk_length};
		mc::arr1d contact_2_upper_corner = {cross_section_1, cross_section_2, total_length};

		_contacts.emplace_back(contact_1_lower_corner, contact_1_upper_corner);
		_contacts.emplace_back(contact_2_lower_corner, contact_2_upper_corner);
		_bulk.set_borders(bulk_lower_corner, bulk_upper_corner);


		_pilot = std::make_shared<mc::free_flight>();
		_scatterer = std::make_shared<mc::scatter>();

		_contacts[0].populate(_beta, 1100, _pilot, _scatterer);
		_contacts[1].populate(_beta, 100, _pilot, _scatterer);
		_bulk.populate(_beta,0, _pilot, _scatterer);

		_number_of_profile_sections = 10;
		_population_profile.first = 0;
		_population_profile.second = std::vector<mc::t_uint>(_number_of_profile_sections, 0);
		_current_profile.first = 0;
		_current_profile.second = std::vector<mc::t_float>(_number_of_profile_sections, 0.);
	};
	void process_command_line_args(int argc=0, char* argv[]=nullptr) // set the output directory and the output file name
	{
		namespace fs = std::experimental::filesystem;

		std::cout << "current path is " << fs::current_path() << std::endl;

		if (argc <= 1)
		{
			if (fs::exists("/Users"))
			{
				_output_directory.assign("/Users/amirhossein/research/test");
			}
			else
			{
				_output_directory.assign("/home/amirhossein/research/test");
			}
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

		if (fs::is_directory(_output_directory.path()))
		{
			if (not fs::is_empty(_output_directory.path()))
			{
				std::cout << "warning: output directory is NOT empty!!!" << std::endl;
				std::cout << "output directory: " << _output_directory.path() << std::endl;
				std::cout << "deleting the existing directory!!!" << std::endl;
				fs::remove_all(_output_directory.path());
				fs::create_directories(_output_directory.path());
				// std::exit(EXIT_FAILURE);
			}
		}
		else
		{
			std::cout << "error: output path is NOT a directory!!!" << std::endl;
			std::cout << "output path: " << _output_directory.path() << std::endl;
			std::exit(EXIT_FAILURE);
		}
	};
	inline mc::t_uint number_of_particles() // returns the number of particles
	{
		mc::t_uint number_of_particles = _bulk.number_of_particles();
		for (auto& contact: _contacts)
		{
			number_of_particles += contact.number_of_particles();
		}
		return number_of_particles;
	};
	inline void step(mc::t_float dt = 0.1) // step the simulation in time
	{
		for (auto& contact: _contacts)
		{
			for (auto&& it=contact.particles().begin(); it!=contact.particles().end(); ++it)
			{
				it->step(dt,_domain);
				_bulk.enlist(it,contact);
			}
		}

		// move bulk to left_contact and right_contact
		for (auto&& it=_bulk.particles().begin(); it!=_bulk.particles().end(); ++it)
		{
			it->step(dt,_domain);
			for (auto& contact: _contacts)
			{
				bool enlisted = contact.enlist(it,_bulk);
				if (enlisted)
				{
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
	inline void write_state(std::fstream& file) // write the state of the simulation in the output file
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
	inline const mc::t_float& time() const // get the mc simulation time
	{
		return _time;
	};
	inline void repopulate_contacts() // repopulate contacts
	{
		_contacts[0].populate(_beta, 1100, _pilot, _scatterer);
		_contacts[1].populate(_beta, 100, _pilot, _scatterer);
	};
	inline void population_profiler(const mc::t_uint& max_history, std::fstream& file, std::fstream& debug_file) // calculate and save the population profile
	{
		if (_population_profile.first < max_history)
		{
			// open the debug file
			if (! debug_file.is_open())
			{
				debug_file.open(_output_directory.path() / "debug.dat", std::ios::out);
				file << std::showpos << std::scientific;
			}


			_population_profile.first += 1;

			mc::t_float length = (_bulk.upper_corner(2)-_bulk.lower_corner(2))/mc::t_float(_population_profile.second.size());
			int i;
			for (const auto& m_particle: _bulk.particles())
			{
				i = std::floor((m_particle.pos(2)-_bulk.lower_corner(2))/length);
				if ((i >= _population_profile.second.size()) or (i<0))
				{
					debug_file << "index out of bound when making population profile\n";
					debug_file << "particle position = " << m_particle.pos(2) << " , bulk limits = [ " << _bulk.lower_corner(2) << " , " << _bulk.upper_corner(2) << " ]" << std::endl;
				}
				else
				{
					_population_profile.second[i] += 1;
				}
			}
		}
		else
		{
			mc::t_float section_volume = _bulk.volume()/mc::t_float(_population_profile.second.size());

			if (! file.is_open())
			{
				file.open(_output_directory.path() / "population_profile.dat", std::ios::out);
				file << std::showpos << std::scientific;

				// store the position of the middle point of each section
				file << time() << " ";
				mc::t_float length = (_bulk.upper_corner(2)-_bulk.lower_corner(2))/mc::t_float(_population_profile.second.size());
				for (int i=0; i<_population_profile.second.size(); ++i)
				{
					file << mc::t_float(i)*length << " ";
				}
				file << std::endl;
			}

			file << time() << " ";
			for (auto& element: _population_profile.second)
			{
				file << mc::t_float(element)/mc::t_float(_population_profile.first)/section_volume << " ";
				element = 0;
			}
			file << std::endl;

			_population_profile.first = 0;
		}
	};
	inline void save_region_current(const mc::t_uint& max_history, std::fstream& current_file, const mc::t_float& time_step) // save the net currents for each region calculated by counting in and out flow of particles in each contact
	{
		if ( _history_of_region_currents < max_history)
		{
			_history_of_region_currents += 1;
		}
		else
		{
			if (! current_file.is_open())
			{
				current_file.open(_output_directory.path() / "region_current.dat", std::ios::out);
				current_file << std::showpos << std::scientific;
			}

			mc::t_float cross_section = (_bulk.upper_corner(0)-_bulk.lower_corner(0)) * (_bulk.upper_corner(1)-_bulk.lower_corner(1));
			mc::t_float elapsed_time = mc::t_float(_history_of_region_currents)*time_step;

			current_file << time() << " ";
			for (auto& contact: _contacts)
			{
				current_file << mc::q0*mc::t_float(contact.particle_flow())/elapsed_time/cross_section;
				current_file << " ";
				contact.reset_particle_flow();
			}
			current_file << std::endl;

			_history_of_region_currents = 0;
		}
	};
	inline void current_profiler(const mc::t_uint& max_history, std::fstream& file) // calculate and save the current profile using the average velocity method
	{
		if (_current_profile.first < max_history)
		{
			_current_profile.first += 1;

			mc::t_float length = (_bulk.upper_corner(2)-_bulk.lower_corner(2))/mc::t_float(_current_profile.second.size());
			int i;
			for (const auto& m_particle: _bulk.particles())
			{
				i = std::floor((m_particle.pos(2)-_bulk.lower_corner(2))/length);
				_current_profile.second[i] += m_particle.velocity(2);
			}
		}
		else
		{
			mc::t_float section_volume = _bulk.volume()/mc::t_float(_current_profile.second.size());

			if (! file.is_open())
			{
				file.open(_output_directory.path() / "current_profile.dat", std::ios::out);
				file << std::showpos << std::scientific;

				// store the position of the middle point of each section
				file << time() << " ";
				for (int i=0; i<_current_profile.second.size(); ++i)
				{
					mc::t_float length = (_bulk.upper_corner(2)-_bulk.lower_corner(2))/mc::t_float(_current_profile.second.size());
					file << mc::t_float(i)*length << " ";
				}
				file << std::endl;
			}
			file << time() << " ";
			for (auto& element: _current_profile.second)
			{
				file << mc::q0*element/mc::t_float(_current_profile.first)/section_volume;
				file << " ";
				element = 0;
			}
			file << std::endl;

			_current_profile.first = 0;
		}
	};

}; // end class monte_carlo

} // end namespace mc

#endif
