// main file defining the monte carlo simulation

#include <stdio.h>
#include <chrono>
#include <ctime>

#include "utility.h"
#include "monte_carlo.h"
#include "scatter.h"

int main(int argc, char *argv[])
{

	// print the start time and start recording the run time
	std::time_t start_time = std::time(nullptr);
	std::cout << "\n***\nstart time:\n" <<  std::asctime(std::localtime(&start_time)) << "***\n\n";

	mc::monte_carlo mc_simulation;
	mc_simulation.process_command_line_args(argc, argv);

	std::fstream population_file;
	std::fstream current_profile_file;
	std::fstream region_current_file;
	std::fstream debug_file;

	mc::t_float time_step = 1.e-14;
	mc::t_uint max_history = 1000;

	// while (mc_simulation.time() < 5.e-10)
	while (true)
	{
		mc_simulation.step(time_step);
		mc_simulation.repopulate_contacts();
		mc_simulation.population_profiler(max_history, population_file, debug_file);
		mc_simulation.save_region_current(max_history, region_current_file, time_step);

		if (int(mc_simulation.time() / time_step) % 1000 == 0)
		{
			std::cout << "simulation time [seconds]: " << mc_simulation.time() << " .... number of particles: " << mc_simulation.number_of_particles() << std::endl;
		}
	}

	// print the end time and the runtime
	std::time_t end_time = std::time(nullptr);
	std::cout << "\nruntime: " << std::difftime(end_time,start_time) << " seconds" << std::endl;
	std::cout << "\n***\nend time:\n" << std::asctime(std::localtime(&end_time)) << "***\n\n";

	return 0;
}
