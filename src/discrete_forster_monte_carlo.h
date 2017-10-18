#ifndef discrete_forster_monte_carlo_h
#define discrete_forster_monte_carlo_h

#include <iostream>
#include <list>
#include <experimental/filesystem>
#include <fstream>
#include <map>
#include <chrono>
#include <cmath>

#include "../rapidxml/rapidxml.hpp"
#include "../rapidxml/rapidxml_utils.hpp"
#include "../rapidxml/rapidxml_print.hpp"

#include "all_particles.h"
#include "utility.h"
#include "region.h"
#include "discrete_forster_region.h"

namespace mc
{

class discrete_forster_monte_carlo
{

public:
	typedef mc::discrete_forster_particle t_particle;
	typedef mc::discrete_forster_free_flight t_ff;
	typedef mc::discrete_forster_scatter t_scatter;
	typedef mc::discrete_forster_region t_region;

private:
	mc::t_float _time;
	mc::t_float _max_hopping_rate;

	std::vector<t_region> _regions;
	std::list<std::shared_ptr<t_scatter>> _all_scat_list;
	std::pair<mc::arr1d, mc::arr1d> _domain;

	mc::t_uint _number_of_contact1_particles, _number_of_contact2_particles; // number of particles in the contacts

	std::experimental::filesystem::directory_entry _output_directory; // this is the address of the output_directory
	std::experimental::filesystem::directory_entry _input_directory; // this is the address of the output_directory
	std::tuple<std::vector<mc::t_uint>, mc::t_uint, mc::t_float, mc::t_float> _population_profile; // this is the population profile of particles through the simulation domain along the z-axis
	mc::t_uint _history_of_region_currents; // this is the number of steps that the net _current in the regions have been recorded

public:

	discrete_forster_monte_carlo() {};
	// get the mc simulation time
	const mc::t_float& time() const
	{
		return _time;
	};
	// increase simulation time by dt
	void increase_time(const mc::t_float& dt)
	{
		_time += dt;
	};
	// get constant reference to the output_directory
	const std::experimental::filesystem::directory_entry& output_directory()
	{
		return _output_directory;
	};
	// get constant reference to the input_directory
	const std::experimental::filesystem::directory_entry& input_directory()
	{
		return _input_directory;
	};
	// get constant reference to the output_directory
	const std::experimental::filesystem::path& output_path()
	{
		return _output_directory.path();
	};
	// get constant reference to the input_directory
	const std::experimental::filesystem::path& input_path()
	{
		return _input_directory.path();
	};
	// returns the number of particles
	mc::t_uint number_of_particles()
	{
		mc::t_uint number_of_particles = 0;
		for (const auto& m_region: _regions)
		{
			number_of_particles += m_region.number_of_particles();
		}
		return number_of_particles;
	};
	//initialize the simulation condition
	void init()
	{
		create_scatterers(input_directory().path(), output_directory().path());
		find_neighbors(_max_hopping_rate, 40.e-9);
		_domain = find_minmax_coordinates();

		mc::t_float x_min = _domain.first[0];
		mc::t_float x_max = _domain.second[0];

		mc::t_float y_min = _domain.first[1];
		mc::t_float y_1   = _domain.first[1]+0.1*(_domain.second[1]-_domain.first[1]);
		mc::t_float y_2   = _domain.first[1]+0.9*(_domain.second[1]-_domain.first[1]);
		mc::t_float y_max = _domain.second[1];

		mc::t_float z_min = _domain.first[2];
		mc::t_float z_max = _domain.second[2];

		mc::arr1d contact_1_lower_corner = {x_min, y_min, z_min};
		mc::arr1d contact_1_upper_corner = {x_max, y_1,   z_max};
		mc::arr1d      bulk_lower_corner = {x_min, y_1,   z_min};
		mc::arr1d      bulk_upper_corner = {x_max, y_2,   z_max};
		mc::arr1d contact_2_lower_corner = {x_min, y_2,   z_min};
		mc::arr1d contact_2_upper_corner = {x_max, y_max, z_max};

		for (int i=0; i<3; i++)
		{
			_regions.push_back(mc::discrete_forster_region());
		}

		_regions[0].set_borders(contact_1_lower_corner, contact_1_upper_corner);
		_regions[1].set_borders(bulk_lower_corner, bulk_upper_corner);
		_regions[2].set_borders(contact_2_lower_corner, contact_2_upper_corner);

		for (auto& m_region: _regions)
		{
			m_region.create_scatterer_vector(_all_scat_list);
			m_region.create_pilot_list();
		}


		_number_of_contact1_particles = 1100;
		_number_of_contact2_particles = 100;
		_regions.front().populate(_number_of_contact1_particles);
		_regions.back().populate(_number_of_contact2_particles);

		t_region& bulk = _regions[1];
		mc::t_uint number_of_profile_sections = 10;
		std::get<0>(_population_profile) = std::vector<mc::t_uint>(number_of_profile_sections, 0);
		std::get<1>(_population_profile) = 0; // history of the population_profiler
		std::get<2>(_population_profile) = (bulk.upper_corner(1)-bulk.lower_corner(1))/mc::t_float(number_of_profile_sections); // length of each segment in population profiler
		std::get<3>(_population_profile) = bulk.volume()/mc::t_float(number_of_profile_sections); // volume of each segment in population profiler

	};
	// read in the coordinate of all the cnt segments or molecules and create the scatterer objects that manage particle hopping between the sites
	void create_scatterers(const std::experimental::filesystem::path& input_path, const std::experimental::filesystem::path& output_path)
	{
		std::cout << "this is the input path: " << input_path << std::endl;
		std::ifstream file;
		std::string line;

		int num = 1;
		int max_num = 3;
		std::string base = "tube";
		std::string extension = ".dat";
		std::string filename = input_path / (base+std::to_string(num)+extension);

		file.open(filename);
		std::regex tube_rgx("tube");

		// loop over files
		while ((file.is_open()) and (num<max_num))
		{
			std::cout << "reading data from file..." << filename << "...\n";

			while(std::getline(file, line, ';'))
			{

				// find the new tube coordinates by finding the keywork "tube"
				if (!std::regex_search(line,tube_rgx))
				{
					try
					{
						std::istringstream iss(line);
						mc::arr1d pos;
						std::string token;
						int i = 0;
						while(std::getline(iss,token,','))
						{
							pos[i] = 1.e-9*std::stod(token);
							i++;
						}
						if (i!= 3)
						{
							throw std::range_error("Could not read compelete set of coordinates!");
						}
						_all_scat_list.push_back(std::make_shared<mc::discrete_forster_scatter>());
						_all_scat_list.back()->set_pos(pos);
					}
					catch(const std::invalid_argument& e)
					{
						std::cout << e.what() << std::endl;
					}
					catch(const std::range_error& e)
					{
						std::cout << e.what() << std::endl;
					}

				}
			}

			file.close();
			num++;
			filename = input_path / (base+std::to_string(num)+extension);
			file.open(filename);

		}

		std::cout << "total number of discrete_forster_scatterer: " << _all_scat_list.size() << std::endl;

	};
	// find the neighbors of each scattering object
	void find_neighbors(const mc::t_float& max_hopping_radius, const mc::t_float& max_search_radius)
	{
		int counter = 0;

		for (auto i = _all_scat_list.begin(); i != _all_scat_list.end(); ++i)
		{
			for (auto j = std::next(i); j!= _all_scat_list.end(); ++j)
			{
				// mc::t_float distance = std::hypot((*i)->pos(0)-(*j)->pos(0), (*i)->pos(1)-(*j)->pos(1), (*i)->pos(2)-(*j)->pos(2));
				mc::t_float distance = std::sqrt(std::pow((*i)->pos(0)-(*j)->pos(0),2)+ std::pow((*i)->pos(1)-(*j)->pos(1),2) + std::pow((*i)->pos(2)-(*j)->pos(2),2));
				if (distance < max_hopping_radius)
				{
					// std::cout << "found new neighboring: distance = " << distance << "\n";
					(*i)->add_neighbor(*j, distance);
					(*j)->add_neighbor(*i, distance);
				}

				// break the search if the scatterers are getting too far apart to speed up the search process
				if (distance > max_hopping_radius)
				{
					break;
				}
			}
			(*i)->sort_neighbors();
			(*i)->make_cumulative_scat_rate();

			counter ++;
			if (counter %1000 == 0)
			{
				std::cout << "number of neighbors: " << (*i)->number_of_neighbors() << "\n";
				std::cout << counter << "\n";
			}

		}

		std::cout << "finished finding neighbors" << std::endl;
	};
	// find minimum of the minimum coordinates of the scattering objects
	std::pair<mc::arr1d, mc::arr1d> find_minmax_coordinates()
	{
		std::pair<mc::arr1d, mc::arr1d> minmax_coordinates;

		for (int i=0; i<3; i++)
		{
			auto it = _all_scat_list.begin();
			minmax_coordinates.first[i] = (*it)->pos(i);
			minmax_coordinates.second[i] = (*it)->pos(i);

			while(it != _all_scat_list.end())
			{
				if (minmax_coordinates.first[i] > (*it)->pos(i))
				{
					minmax_coordinates.first[i] = (*it)->pos(i);
				}

				if (minmax_coordinates.second[i] < (*it)->pos(i))
				{
					minmax_coordinates.second[i] = (*it)->pos(i);
				}
				it++;
			}
		}

		for (int i=0; i<3; i++)
		{
			std::cout << "min coordinate: " << minmax_coordinates.first[i] << "  ,  max coordinate: " << minmax_coordinates.second[i] << std::endl;
		}

		return minmax_coordinates;
	}
	// step the simulation in time
	void step(mc::t_float dt)
	{

		for (auto&& it=_regions[0].particles().begin(); it!=_regions[0].particles().end(); ++it)
		{
			(*it)->step(dt,_domain);
			_regions[1].enlist(it,_regions[0]);
		}


		for (auto&& it=_regions[2].particles().begin(); it!=_regions[2].particles().end(); ++it)
		{
			(*it)->step(dt,_domain);
			_regions[1].enlist(it,_regions[2]);
		}

		for (auto&& it=_regions[1].particles().begin(); it!=_regions[1].particles().end(); ++it)
		{
			(*it)->step(dt,_domain);
			if (not _regions[2].enlist(it,_regions[1]))
			{
				_regions[0].enlist(it,_regions[1]);
			}
		}

		// dump the new particles into the particles list in each region
		for (auto&& m_region: _regions)
		{
			m_region.dump_new_particles();
		}

		increase_time(dt);

	};
	// repopulate contacts
	void repopulate_contacts()
	{
		_regions.front().populate(_number_of_contact1_particles);
		_regions.back().populate(_number_of_contact2_particles);
	};
	// calculate and save the population profile
	void population_profiler(const mc::t_uint& max_history, std::fstream& file, std::fstream& debug_file)
	{

		if (std::get<1>(_population_profile) < max_history)
		{
			// open the debug file
			if (! debug_file.is_open())
			{
				debug_file.open(_output_directory.path() / "debug.dat", std::ios::out);
				file << std::showpos << std::scientific;
			}

			t_region& bulk = _regions[1];

			std::get<1>(_population_profile) += 1;

			int i;
			for (const auto& m_particle_ptr: bulk.particles())
			{
				i = std::floor((m_particle_ptr->pos(1) - bulk.lower_corner(1))/std::get<2>(_population_profile));
				if ((i >= std::get<0>(_population_profile).size()) or (i<0))
				{
					debug_file << "index out of bound when making population profile\n";
					debug_file << "particle position = " << m_particle_ptr->pos(1) << " , bulk limits = [ " << bulk.lower_corner(1) << " , " << bulk.upper_corner(1) << " ]" << std::endl;
				}
				else
				{
					(std::get<0>(_population_profile))[i] += 1;
				}
			}
		}
		else
		{
			if (not file.is_open())
			{
				file.open(_output_directory.path() / "population_profile.dat", std::ios::out);
				file << std::showpos << std::scientific;

				// store the position of the middle point of each section
				file << time() << " ";
				for (int i=0; i < std::get<0>(_population_profile).size(); ++i)
				{
					file << mc::t_float(i)*std::get<2>(_population_profile) << " ";
				}
				file << std::endl;
			}

			file << time() << " ";
			for (auto& element: std::get<0>(_population_profile))
			{
				file << mc::t_float(element)/mc::t_float(std::get<1>(_population_profile))/std::get<3>(_population_profile) << " ";
				element = 0;
			}
			file << std::endl;

			std::get<1>(_population_profile) = 0;
		}
	};
	// save the net currents for each region calculated by counting in and out flow of particles in each contact
	void save_region_current(const mc::t_uint& max_history, std::fstream& current_file, const mc::t_float& time_step)
	{
		if ( _history_of_region_currents < max_history)
		{
			_history_of_region_currents += 1;
		}
		else
		{
			if (! current_file.is_open())
			{
				current_file.open(_output_directory.path() / "region_current.dat", std::ios::out);
				current_file << std::showpos << std::scientific;
			}

			t_region& bulk = _regions[1];
			mc::t_float cross_section = (bulk.upper_corner(0)-bulk.lower_corner(0)) * (bulk.upper_corner(2)-bulk.lower_corner(2));
			mc::t_float elapsed_time = mc::t_float(_history_of_region_currents)*time_step;

			current_file << time() << " ";
			for (auto&& m_region: _regions)
			{
				current_file << mc::t_float(m_region.particle_flow())/elapsed_time/cross_section;
				current_file << " ";
				m_region.reset_particle_flow();
			}
			current_file << std::endl;

			_history_of_region_currents = 0;
		}
	};
	// set the output directory and the output file name
	void process_command_line_args(int argc, char* argv[])
	{
		namespace fs = std::experimental::filesystem;


		// first find the xml input file and open it

		fs::directory_entry xml_file;

		std::cout << "current path is " << fs::current_path() << std::endl;

		if (argc <= 1)
		{
			xml_file.assign("input.xml");
		}
		else
		{
			xml_file.assign(argv[1]);
		}

		if(fs::exists(xml_file))
		{
			std::cout << "input xml file found: " << xml_file.path() << std::endl;
		}
		else
		{
			std::cout << "input xml file NOT found: " << xml_file.path() << std::endl;
			std::exit(1);
		}

		if (!fs::is_regular_file(xml_file))
		{
			std::cout << "input xml file NOT found: " << xml_file.path() << std::endl;
			std::exit(1);
		}
		std::cout << std::endl;

		rapidxml::file<> xmlFile(xml_file.path().c_str()); //open file
		rapidxml::xml_document<> doc; //create xml object
		doc.parse<0>(xmlFile.data()); //parse contents of file
		rapidxml::xml_node<>* curr_node = doc.first_node(); //gets the node "Document" or the root nodes

		// set the output_directory
		{
			curr_node = curr_node->first_node("output_directory");
			std::string attr = curr_node->first_attribute("type")->value();
			std::string path = mc::trim(curr_node->value());
			if (attr == "absolute")
			{
				std::cout << "absolute directory format used!\n";
			}

			_output_directory.assign(path);
			std::cout << "output_directory: " << _output_directory.path() << std::endl;

			if (not fs::exists(_output_directory.path()))
			{
				std::cout << "warning: output directory does NOT exist!!!" << std::endl;
				std::cout << "output directory: " << _output_directory.path() << std::endl;
				fs::create_directories(_output_directory.path());
			}

			if (fs::is_directory(_output_directory.path()))
			{
				if (not fs::is_empty(_output_directory.path()))
				{
					std::cout << "warning: output directory is NOT empty!!!" << std::endl;
					std::cout << "output directory: " << _output_directory.path() << std::endl;
					std::cout << "deleting the existing directory!!!" << std::endl;
					fs::remove_all(_output_directory.path());
					fs::create_directories(_output_directory.path());
				}
			}
			else
			{
				std::cout << "error: output path is NOT a directory!!!" << std::endl;
				std::cout << "output path: " << _output_directory.path() << std::endl;
				std::exit(EXIT_FAILURE);
			}
		}

		// set input directory
		{
			auto next_node = curr_node->next_sibling("input_directory");
			if (next_node == 0)
			{
				next_node = curr_node->previous_sibling("input_directory");
				if (next_node == 0)
				{
					std::cout << "input_directory not found in XML file!!!" << std::endl;
					std::exit(1);
				}
			}
			curr_node = next_node;

			std::string attr = curr_node->first_attribute("type")->value();
			std::string path = mc::trim(curr_node->value());
			if (attr == "absolute")
			{
				std::cout << "absolute directory format used!\n";
			}
			_input_directory.assign(path);

			if (not fs::exists(_input_directory.path()))
			{
				std::cout << "\n***\nwarning: input directory does NOT exist!!!\n"
									<< "input directory: " << _input_directory.path() << "\n***\n\n";
				std::exit(1);
			}

			if (fs::is_directory(_input_directory.path()))
			{
				if (fs::is_empty(_input_directory.path()))
				{
					std::cout << "warning: input directory is empty!!!" << std::endl;
					std::cout << "input directory: " << _input_directory.path() << std::endl;
				}
			}
			else
			{
				std::cout << "\n***\nerror: input path is NOT a directory!!!\n"
									<< "input path: " << _input_directory.path() << std::endl;
				std::exit(EXIT_FAILURE);
			}
		}

		// std::exit(0);

	};

}; // end class discrete_forster_monte_carlo

} // end namespace mc

#endif // discrete_forster_monte_carlo_h
