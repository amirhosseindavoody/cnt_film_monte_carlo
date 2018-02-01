// main file defining the monte carlo simulation

#include <stdio.h>
#include <chrono>
#include <ctime>

#include "./helper/utility.h"
#include "./discrete_forster/discrete_forster_monte_carlo.h"
#include "./lib/json.hpp"

#include "./exciton_transfer/cnt.h"

int main(int argc, char *argv[])
{

	// print the start time and start recording the run time
	std::time_t start_time = std::time(nullptr);
	std::cout << "\n***\nstart time:\n" <<  std::asctime(std::localtime(&start_time)) << "***\n\n";


	using json = nlohmann::json;

	std::string filename;
	if (argc <= 1){
		filename = "input.json";
	} else {
		filename = argv[1];
	}

	// read a JSON file
	std::ifstream input_file(filename.c_str());
	json j;
	input_file >> j;

	// get the json part related to exciton mc simulation
	if (j.count("exciton monte carlo")==0){
		throw std::invalid_argument("json input file does not contain \"exciton monte carlo\"");
	}
	nlohmann::json json_mc = j["exciton monte carlo"];

	// if exciton transfer type is davoody get cnt json information and add it to json_mc
	if (std::string(j["exciton monte carlo"]["rate type"]) == "davoody"){
		json_mc["cnts"] = j["cnts"];
	}

	std::cout << std::setw(4) << json_mc << std::endl;

	throw std::out_of_range("check to see if json_mc is printed correctly before moving on");

	//***********************************************************************************************
	// create cnts and calculate their exciton states
	//***********************************************************************************************
	
	std::vector<cnt> cnts;
	// get the parent directory for cnts
	std::string parent_directory = j["cnts"]["directory"];
	j["cnts"].erase("directory");

	// create excitons and calculate exciton dispersions
	cnts.reserve(j["cnts"].size()); // this is reservation of space is crucial to ensure we do not move cnts, since the move constructor is not implemented yet
	for (const auto& j_cnt: j["cnts"])
	{
		cnts.emplace_back(cnt(j_cnt,parent_directory));
		cnts.back().calculate_exciton_dispersion();
	};

	//***********************************************************************************************
	// create monte carlo object and run the MC simulation
	//***********************************************************************************************

	if (j.count("exciton monte carlo")==0){
		throw std::invalid_argument("input.json should contain \"exciton monte carlo\" property.");
	}
	
	mc::discrete_forster_monte_carlo sim(json_mc);
	// sim.initialize_scattering_table(cnts[0],cnts[0]);


	// initialize and run simulation for the exciton hopping
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

	// print the end time and the runtime
	std::time_t end_time = std::time(nullptr);
	std::cout << "\nruntime: " << std::difftime(end_time,start_time) << " seconds" << std::endl;
	std::cout << "\n***\nend time:\n" << std::asctime(std::localtime(&end_time)) << "***\n\n";

	return 0;
}
