#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <regex>
#include <list>
#include <vector>
#include <memory>

#include "file_management.h"

using namespace std;

file_management:: file_management()
{
}

file_management:: ~file_management()
{
}

// this function sets the input directory and checks if it exists and accessible.
void file_management::set_input_directory(const char *input)
{
	input_directory = input;

	// check if folder or file exists.
	struct stat buf;
	if (stat(input_directory.c_str(),&buf) != 0)
	{
		cout << "Error: incorrect result directory path!!!" << endl;;
		exit(EXIT_FAILURE);
	}

	// check if there is '/' at the end of the input directory
	char last_char = input_directory.at(input_directory.size()-1);
	if (last_char != '/')
	{
		input_directory.push_back('/');
	}
	cout << "input directory: " << input_directory << endl;

}

// this function gives back the input directory
string file_management::get_input_directory()
{
	return input_directory;
}

// This function creates a list of files that contain information of cnt geometry and location from the bullet physics simulation.
void file_management::build_cnt_file_list()
{
	file_list = make_unique<list<string>>();
	DIR *resDir;
	struct dirent *ent = nullptr;

	regex rgx("CNT_Num_\\d+\\.csv"); //files we want to look through

	//Check if folder can be opened - should work due to above checks
	if ((resDir = opendir(input_directory.c_str())) != nullptr)
	{
		//throw away first two results as they are . and ..
		readdir(resDir);
		readdir(resDir);

		//iterate over all of the real files
		while ((ent = readdir(resDir)) != nullptr)
		{
			smatch matches; //match_results for string objects
			string tmps =  string(ent->d_name);
			regex_search(tmps, matches, rgx);
			if (!matches.empty())
			{
				file_list->push_back(ent->d_name);
			}
		}
		closedir(resDir); //deletes pointer

	}
	else
	{
		cout << "Could not open results directory!!!" << endl;
		exit(EXIT_FAILURE);
	}
	delete ent;
	ent = nullptr;
}


// void file_management::parse_xml(string input_path)
// {

// }