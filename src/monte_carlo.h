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
	mc::t_uint _num_particles; // number of particles
	mc::t_float _time; // total simulation time that has elapsed
	mc::t_float _temperature; // temperature of the simulation
	mc::t_float _beta; // 1/kB*T is the inverse of the thermal energy
	mc::arr1d _volume; // the physical extent of the simulation domain

	// defining various domains in the simulation
	mc::region _bulk;
	std::vector<mc::region> _contacts;

	// output directory
	std::experimental::filesystem::directory_entry _output_directory; // this is the address of the output_directory

public:

	monte_carlo(unsigned long int num_particles = 100); // constructor
	void process_command_line_args(int argc=0, char* argv[]=nullptr); // set the output directory and the output file name
	mc::t_uint num_particles(); // returns the number of particles
	void step(mc::t_float dt = 0.1); // step the simulation in time
	void write_state(std::fstream& file); // write the state of the simulation in the output file
	inline const mc::t_float& time() const // get the mc simulation time
	{
		return _time;
	};

}; // end class monte_carlo

} // end namespace mc

#endif
