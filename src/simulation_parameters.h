#ifndef simulation_parameters_h
#define simulation_parameters_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <regex>
#include <list>
#include <vector>
#include <memory>

using namespace std;

class simulation_parameters
{
	public:

		simulation_parameters();

		double segment_length; //[Angstrom] this is the length of the cnt segments that are used in the monte carlo simulation
		double region_length_min; //[Angstrom]
		int number_of_steps; //the total number of time steps in the monte carlo simulation
		double maximum_distance; //[Angstroms]  the maximum segment separation to be considered for exciton transfer
		int number_of_excitons; //the number of excitons in the simulation
		bool auto_complete;
		double threshold; //this is the criteria for ending the monte carlo simulation: difference/max_difference < threshold
		int num_to_finish;// number of differences/maxDiff  that must be below thresh in a row to finish
		int num_to_check; //number of time steps to check completion
		double tfac; //percent free flight times above delta T
	private:
		
		
		
};


#endif // simulation_parameters_h