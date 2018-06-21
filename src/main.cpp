// main file defining the monte carlo simulation

#include <stdio.h>
#include <chrono>
#include <ctime>

#include <experimental/filesystem>

#include "../lib/json.hpp"
#include "./discrete_forster/monte_carlo.h"
#include "./helper/utility.h"

#include "./exciton_transfer/cnt.h"
#include "./helper/prepare_directory.hpp"

int main(int argc, char *argv[])
{

	// print the start time and start recording the run time
	std::time_t start_time = std::time(nullptr);
	std::cout << "\n***\nstart time:\n" <<  std::asctime(std::localtime(&start_time)) << "***\n\n";

	// get the input JSON filename
	std::string filename;
	if (argc <= 1){
		filename = "input.json";
	} else {
		filename = argv[1];
	}

	// read the input JSON file
	std::ifstream input_file(filename.c_str());
	nlohmann::json j;
	input_file >> j;

	// get the json part related to exciton mc simulation
	if (j.count("exciton monte carlo")==0){
		throw std::invalid_argument("json input file does not contain \"exciton monte carlo\"");
	}
	nlohmann::json json_mc = j["exciton monte carlo"];

	// if exciton transfer type is davoody get cnt json information and add it to json_mc
	if (j["exciton monte carlo"]["rate type"].get<std::string>() == "davoody"){
		json_mc["cnts"] = j["cnts"];
	}

	//***********************************************************************************************
	// create monte carlo object and run the MC simulation
	//***********************************************************************************************

	if (j.count("exciton monte carlo")==0){
		throw std::invalid_argument("input.json should contain \"exciton monte carlo\" property.");
	}
	
	mc::monte_carlo sim(json_mc);

	// initialize and run simulation for the exciton hopping
	sim.init();
	sim.save_json_properties();

	double time_step = 1.e-14;

  std::cout << "running Monte Carlo:\n";

	while (true)
	{
		
    // std::cout << "simulation time [seconds]: " << std::scientific << sim.time()
    //           << " .... "
    //           << "number of particles: " << sim.number_of_particles()
    //           << std::endl;

    sim.step(time_step);
    sim.save_metrics();

		sim.repopulate_contacts();


    std::cout << "simulation time [seconds]: " << std::scientific << sim.time()
              << " .... "
              << "number of particles: " << sim.number_of_particles()
              << "\r" << std::flush;
	}

	// print the end time and the runtime
	std::time_t end_time = std::time(nullptr);
	std::cout << "\nruntime: " << std::difftime(end_time,start_time) << " seconds" << std::endl;
	std::cout << "\n***\nend time:\n" << std::asctime(std::localtime(&end_time)) << "***\n\n";

	return 0;
}
