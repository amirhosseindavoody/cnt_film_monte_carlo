#ifndef file_management_h
#define file_management_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <regex>
#include <list>
#include <vector>
#include <memory>

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "simulation_parameters.h"

using namespace std;

class file_management
{
	public:
		file_management();
		~file_management();
		void set_input_directory(const char *input); // this function sets the input directory and checks if it exists and accessible.
		string get_input_directory(); // this function gives back the input directory
		void build_cnt_file_list();  // This function creates a list of files that contain information of cnt geometry and location from the bullet physics simulation.
		simulation_parameters parse_xml(); // this function parses the input xml file and puts the simulation parameters in the simulation parameter object "sim".
		void change_working_directory(const char *path);




		string output_directory; // this is the address of the output folder that the simulation results would be saved into.
		unique_ptr<list<string>>  file_list; // this list contains a list of the files that contain information of the geometry and location of CNTs.

	private:
		string input_directory; // this is the address of the input folder that contains information of location and geometry of CNTs.

		
};


#endif // file_management_h