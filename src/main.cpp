// main file defining the monte carlo simulation

#include <stdio.h>
#include <chrono>

#include "utility.h"
#include "monte_carlo.h"
#include "scatter.h"

int main(int argc, char *argv[])
{

	auto start_time = std::chrono::high_resolution_clock::now();

	mc::monte_carlo mc_simulation(2);
	mc_simulation.process_command_line_args(argc, argv);

	float time_step = 1.e-14;

	std::fstream file;

	while (mc_simulation.time() < 1.e-10)
	{
		// std::cout << "\nsimulation time [seconds]: " << mc_simulation.time() << std::endl << std::endl;

		mc_simulation.step(time_step);

		if (int(mc_simulation.time() / time_step) % 1000 == 0)
		{
			std::cout << "simulation time [seconds]: " << mc_simulation.time() << std::endl;
			// mc_simulation.update_particle_list();
			// mc_simulation.write_state(file);
		}
	}

	// mc_simulation.write_state(file);

	auto end_time = std::chrono::high_resolution_clock::now();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
	std::cout << std::scientific << "runtime is: " << float(elapsed_time.count())/1.e3 << " [seconds]" << std::endl;

	return 0;
}
