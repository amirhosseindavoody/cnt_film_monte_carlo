// main file defining the monte carlo simulation

#include <stdio.h>
#include <chrono>
#include <ctime>

#include "./helper/utility.h"
#include "./discrete_forster/discrete_forster_monte_carlo.h"

int main(int argc, char *argv[])
{

	// print the start time and start recording the run time
	std::time_t start_time = std::time(nullptr);
	std::cout << "\n***\nstart time:\n" <<  std::asctime(std::localtime(&start_time)) << "***\n\n";


	// initialize and run simulation for the exciton hopping
	mc::discrete_forster_monte_carlo sim;
	sim.process_command_line_args(argc, argv);
	sim.init();

	std::fstream population_file;
	std::fstream current_profile_file;
	std::fstream region_current_file;
	std::fstream debug_file;

	mc::t_float time_step = 1.e-14;
	mc::t_uint max_history = 100;

	while (true)
	{
		sim.step(time_step);
		sim.repopulate_contacts();
		sim.population_profiler(max_history, population_file, debug_file);
		sim.save_region_current(max_history, region_current_file, time_step);

		if (int(sim.time() / time_step) % 10 == 0)
		{
			std::cout << "simulation time [seconds]: " << std::scientific << sim.time() << " .... number of particles: " << sim.number_of_particles() <<"\r" << std::flush;
		}
	}

	// // // initialize and run simulation for the particle diffusion model
	// mc::monte_carlo sim;
	// sim.process_command_line_args(argc, argv);
	//
	// std::fstream population_file;
	// std::fstream current_profile_file;
	// std::fstream region_current_file;
	// std::fstream debug_file;
	//
	// mc::t_float time_step = 1.e-14;
	// mc::t_uint max_history = 1000;
	//
	// // while (sim.time() < 5.e-10)
	// while (true)
	// {
	// 	sim.step(time_step);
	// 	sim.repopulate_contacts();
	// 	sim.population_profiler(max_history, population_file, debug_file);
	// 	sim.save_region_current(max_history, region_current_file, time_step);
	//
	// 	if (int(sim.time() / time_step) % 1000 == 0)
	// 	{
	// 		std::cout << "simulation time [seconds]: " << sim.time() << " .... number of particles: " << sim.number_of_particles() << std::endl;
	// 	}
	// }

	// print the end time and the runtime
	std::time_t end_time = std::time(nullptr);
	std::cout << "\nruntime: " << std::difftime(end_time,start_time) << " seconds" << std::endl;
	std::cout << "\n***\nend time:\n" << std::asctime(std::localtime(&end_time)) << "***\n\n";

	return 0;
}
