#ifndef monte_carlo_h
#define monte_carlo_h

#include <iostream>
#include <list>
#include <experimental/filesystem>
#include <fstream>
#include <map>
#include <chrono>

#include "all_particles.h"

#include "utility.h"
#include "region.h"

namespace mc
{

class monte_carlo
{
private:

	// simulation parameters
	mc::t_uint _num_particles; // number of particles
	mc::t_float _temperature; // temperature of the simulation
	mc::t_float _beta; // 1/kB*T is the inverse of the thermal energy
	mc::t_float _time; // total simulation time that has elapsed
	std::pair<mc::arr1d, mc::arr1d> _domain; // the physical extent of the simulation domain

	std::shared_ptr<mc::free_flight> _pilot; // pilot for the free flight step for particles in the simulation
	std::shared_ptr<mc::scatter> _scatterer; // scatterer for all particles in the simulation

	// defining various domains in the simulation
	mc::region _bulk;
	std::vector<mc::region> _contacts;

	// output files properties
	std::experimental::filesystem::directory_entry _output_directory; // this is the address of the output_directory
	std::experimental::filesystem::directory_entry _input_directory; // this is the address of the output_directory
	mc::t_uint _number_of_profile_sections; // number of profile sections used for profiling
	std::pair<mc::t_uint, std::vector<mc::t_uint>> _population_profile; // this is the population profile of particles through the simulation domain along the z-axis
	std::pair<mc::t_uint, std::vector<mc::t_float>> _current_profile; // this is the average current profile along the simulation domain
	mc::t_uint _history_of_region_currents; // this is the number of steps that the net _current in the regions have been recorded

protected:

public:

	monte_carlo(); // constructor
	void process_command_line_args(int argc=0, char* argv[]=nullptr); // set the output directory and the output file name
	mc::t_uint number_of_particles() // returns the number of particles
	{
		mc::t_uint number_of_particles = _bulk.number_of_particles();
		for (auto& contact: _contacts)
		{
			number_of_particles += contact.number_of_particles();
		}
		return number_of_particles;
	};
	void step(mc::t_float dt) // step the simulation in time
	{
		for (auto& contact: _contacts)
		{
			for (auto&& it=contact.particles().begin(); it!=contact.particles().end(); ++it)
			{
				(*it)->step(dt,_domain);
				_bulk.enlist(it,contact);
			}
		}

		// move bulk to left_contact and right_contact
		for (auto&& it=_bulk.particles().begin(); it!=_bulk.particles().end(); ++it)
		{
			(*it)->step(dt,_domain);
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

		increase_time(dt);

	};
	inline void write_state(std::fstream& file) // write the state of the simulation in the output file
	{
		if (! file.is_open())
		{
			file.open(_output_directory.path() / "output.dat", std::ios::out);
			file << std::showpos << std::scientific;
		}
	};
	inline const mc::t_float& time() const // get the mc simulation time
	{
		return _time;
	};
	inline void increase_time(const mc::t_float& dt) // increase simulation time by dt
	{
		_time += dt;
	}
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
			for (const auto& m_particle_ptr: _bulk.particles())
			{
				i = std::floor((m_particle_ptr->pos(2)-_bulk.lower_corner(2))/length);
				if ((i >= _population_profile.second.size()) or (i<0))
				{
					debug_file << "index out of bound when making population profile\n";
					debug_file << "particle position = " << m_particle_ptr->pos(2) << " , bulk limits = [ " << _bulk.lower_corner(2) << " , " << _bulk.upper_corner(2) << " ]" << std::endl;
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
			for (const auto& m_particle_ptr: _bulk.particles())
			{
				i = std::floor((m_particle_ptr->pos(2)-_bulk.lower_corner(2))/length);
				_current_profile.second[i] += m_particle_ptr->velocity(2);
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
	inline void set_temperature(const mc::t_float& value)
	{
		_temperature = value;
		_beta = 1./_temperature/mc::kB;
		// std::cout << "temperature set to " << _temperature << " Kelvins\n";
	};
	const std::experimental::filesystem::directory_entry& output_directory() // get constant reference to the output_directory
	{
		return _output_directory;
	};
	const std::experimental::filesystem::directory_entry& input_directory() // get constant reference to the input_directory
	{
		return _input_directory;
	};
}; // end class monte_carlo

} // end namespace mc

#endif
