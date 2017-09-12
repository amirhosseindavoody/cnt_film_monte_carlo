#ifndef monte_carlo_h
#define monte_carlo_h

#include <iostream>
#include <list>
#include <experimental/filesystem>
#include <fstream>

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
	std::vector<mc::t_uint> _population_profile; // this is the population profile of particles through the simulation domain along the z-axis
	mc::t_uint _history_of_population_profiler; // this is the number of steps that the _population_profile have been recorded for considering the running average
	std::vector<mc::t_float> _current_profile; // this is the average current profile along the simulation domain

public:

	monte_carlo(unsigned long int num_particles = 100); // constructor
	void process_command_line_args(int argc=0, char* argv[]=nullptr); // set the output directory and the output file name
	inline mc::t_uint number_of_particles() // returns the number of particles
	{
		mc::t_uint number_of_particles = _bulk.number_of_particles();
		for (auto& contact: _contacts)
		{
			number_of_particles += contact.number_of_particles();
		}
		return number_of_particles;
	};
	void step(mc::t_float dt = 0.1); // step the simulation in time
	void write_state(std::fstream& file); // write the state of the simulation in the output file
	inline const mc::t_float& time() const // get the mc simulation time
	{
		return _time;
	};
	void update_profile(const mc::t_uint& max_history, std::fstream& population_file, std::fstream& current_file); // calculate and save the population profile

	inline void repopulate_contacts() // repopulate contacts
	{
		_contacts[0].populate(_beta, 1100, _pilot, _scatterer);
		_contacts[1].populate(_beta, 100, _pilot, _scatterer);
	};

}; // end class monte_carlo

} // end namespace mc

#endif
