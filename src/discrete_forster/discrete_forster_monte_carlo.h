#ifndef discrete_forster_monte_carlo_h
#define discrete_forster_monte_carlo_h

#include <iostream>
#include <list>
#include <experimental/filesystem>
#include <fstream>
#include <map>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <algorithm>

#include "../helper/utility.h"
#include "../helper/prepare_directory.hpp"
#include "../helper/constants.h"
#include "../helper/progress.hpp"

#include "../lib/rapidxml/rapidxml.hpp"
#include "../lib/rapidxml/rapidxml_utils.hpp"
#include "../lib/rapidxml/rapidxml_print.hpp"
#include "../lib/json.hpp"


#include "./discrete_forster_region.h"
#include "../exciton_transfer/cnt.h"
#include "../exciton_transfer/exciton_transfer.h"

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

	struct population_profile
	{
		std::vector<unsigned> profile; // number of particles in each segment
		unsigned number_of_sections; // number of sections in the profiler
		unsigned history; // history of the population_profiler
		double dL; // length of each segment in population profiler
		double dV; // volume of each segment in population profiler
		unsigned dim; // this is the dimension along which the population profiler is working
	};

	double _time;
	double _max_hopping_radius;

	nlohmann::json _json_prop; // input properties of the whole mc simulation in json format
	nlohmann::json _json_scat; // input properties of the scattering mechanism in json format

	enum scattering_type {davoody, forster, wong}; // enumerated type to specify type of transfer rate
	scattering_type _scat_t;

	std::array<t_region,3> _regions; // simulation regions: contacts and bulk. Since t_region data type has unique_ptr as data members, we can't use vectors and dynamically push_back new regions. If we want to do that we need to implement the move and copy constructors explicitly.
	std::list<std::shared_ptr<t_scatter>> _all_scat_list; // holds the list of all scattering sites containing the position and orientation
	std::pair<arma::vec, arma::vec> _domain;

	unsigned _number_of_contact1_particles, _number_of_contact2_particles; // number of particles in the contacts

	std::experimental::filesystem::directory_entry _output_directory; // this is the address of the output_directory
	std::experimental::filesystem::directory_entry _input_directory; // this is the address of the input_directory
	unsigned _history_of_region_currents; // this is the number of steps that the net _current in the regions have been recorded
	population_profile _population_probe; // this is the population profile of particles through the simulation domain along the z-axis

	// struct to bundle information related to different scattering mechanisms
	struct scattering_struct
	{
		// default constructor
		scattering_struct() {};

		// constructor
		scattering_struct(const arma::field<arma::cube>& m_rate, const arma::vec& m_theta, const arma::vec& m_z_shift, const arma::vec& m_axis_shift_1, const arma::vec& m_axis_shift_2)
		{
			rate = m_rate;
			theta = m_theta;
			z_shift = m_z_shift;
			axis_shift_1 = m_axis_shift_1;
			axis_shift_2 = m_axis_shift_2;
		};
		arma::field<arma::cube> rate;
		arma::vec theta;
		arma::vec z_shift;
		arma::vec axis_shift_1;
		arma::vec axis_shift_2;
	};
	scattering_struct _scat_table;

public:
	// default constructor
	discrete_forster_monte_carlo() {};
	
	// constructure with json input file
	discrete_forster_monte_carlo(const nlohmann::json& j)
	{
		// store the json properties for use in other methods
		_json_prop = j;

		// set the output directory
		std::string directory_path = j["output directory"];
		bool keep_old_data = true;
		if (j.count("keep old results")==1){
			keep_old_data = j["keep old results"];
		}
		_output_directory = prepare_directory(directory_path,keep_old_data);

		// set the input directory for mesh information
		directory_path = j["mesh input directory"];
		_input_directory = check_directory(directory_path,false);

		// set json info about scattering mechanism
		_json_scat = j;
		_json_scat.erase("mesh input directory");
		_json_scat.erase("output directory");
		_json_scat.erase("keep old results");

		// set maximum hopping radius
		_max_hopping_radius = _json_scat["max hopping radius"];
		std::cout << "maximum hopping radius: " << _max_hopping_radius*1.e9 << " [nm]\n";

		// specify type of transfer rate that is going to be used
		if (j.count("rate type") == 0){
			throw std::invalid_argument("\"rate type\" must be specifieced and must be one of the following: \"davoody\", \"forster\", \"wong\"");
		}
		if (std::string(j["rate type"]) == "davoody"){
			_scat_t = davoody;
		}
		else if (std::string(j["rate type"]) == "forster"){
			_scat_t = forster;
		}
		else if (std::string(j["rate type"]) == "wong"){
			_scat_t = wong;
		}
		else{
			throw std::invalid_argument("\"rate type\" must be specifieced and must be one of the following: \"davoody\", \"forster\", \"wong\"");
		}

	};

	// get the mc simulation time
	const double& time() const
	{
		return _time;
	};

	// increase simulation time by dt
	void increase_time(const double& dt)
	{
		_time += dt;
	};

	// get constant reference to the output_directory
	const std::experimental::filesystem::directory_entry& output_directory() const
	{
		return _output_directory;
	};

	// get constant reference to the input_directory
	const std::experimental::filesystem::directory_entry& input_directory() const
	{
		return _input_directory;
	};

	// get constant reference to the output_directory
	const std::experimental::filesystem::path& output_path() const
	{
		return _output_directory.path();
	};

	// get constant reference to the input_directory
	const std::experimental::filesystem::path& input_path() const
	{
		return _input_directory.path();
	};

	// returns the number of particles
	unsigned number_of_particles() const
	{
		unsigned number_of_particles = 0;
		for (const auto& m_region: _regions)
		{
			number_of_particles += m_region.number_of_particles();
		}
		return number_of_particles;
	};

	//initialize the simulation condition
	void init()
	{
		// read the information of cnt mesh and create a list of all discerete_forster_scatter objects that manage all scattering objects
		create_scatterers_with_orientation(input_directory().path(), output_directory().path());

		_domain = find_minmax_coordinates();

		double x_min = _domain.first[0];
		double x_max = _domain.second[0];

		double y_min = _domain.first[1];
		double y_1   = _domain.first[1]+0.1*(_domain.second[1]-_domain.first[1]);
		double y_2   = _domain.first[1]+0.9*(_domain.second[1]-_domain.first[1]);
		double y_max = _domain.second[1];

		double z_min = _domain.first[2];
		double z_max = _domain.second[2];

		arma::vec contact_1_lower_corner = {x_min, y_min, z_min};
		arma::vec contact_1_upper_corner = {x_max, y_1,   z_max};
		arma::vec      bulk_lower_corner = {x_min, y_1,   z_min};
		arma::vec      bulk_upper_corner = {x_max, y_2,   z_max};
		arma::vec contact_2_lower_corner = {x_min, y_2,   z_min};
		arma::vec contact_2_upper_corner = {x_max, y_max, z_max};

		double max_search_radius = 40.e-9;
		std::cout << "maximum hopping radius: " << _max_hopping_radius*1.e9 << " [nm]\n";
		find_neighbors(_max_hopping_radius, max_search_radius);

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

		_population_probe.number_of_sections = 10;
		_population_probe.dim = 1;
		_population_probe.profile = std::vector<unsigned>(_population_probe.number_of_sections, 0);
		_population_probe.history = 0;
		_population_probe.dL = (bulk.upper_corner(_population_probe.dim)-bulk.lower_corner(_population_probe.dim))/double(_population_probe.number_of_sections);
		_population_probe.dV = bulk.volume()/(_population_probe.number_of_sections);

	};

	// find the neighbors of each scattering object
	void find_neighbors(const double& max_hopping_radius, const double& max_search_radius);

	// find minimum of the minimum coordinates of the scattering objects
	std::pair<arma::vec, arma::vec> find_minmax_coordinates()
	{
		std::pair<arma::vec, arma::vec> minmax_coordinates;

		for (int i=0; i<3; i++)
		{
			auto it = _all_scat_list.begin();
			(minmax_coordinates.first)(i) = (*it)->pos(i);
			(minmax_coordinates.second)(i) = (*it)->pos(i);

			while(it != _all_scat_list.end())
			{
				if ((minmax_coordinates.first)(i) > (*it)->pos(i))
				{
					(minmax_coordinates.first)(i) = (*it)->pos(i);
				}

				if ((minmax_coordinates.second)(i) < (*it)->pos(i))
				{
					(minmax_coordinates.second)(i) = (*it)->pos(i);
				}
				it++;
			}
		}

		std::ios::fmtflags f(std::cout.flags());

		std::cout << "\n simulation domain:\n";
		std::cout << "    x (" << std::fixed << std::showpos << minmax_coordinates.first[0]*1e9 << " , " << (minmax_coordinates.second)(0)*1e9 << ") [nm]\n";
		std::cout << "    y (" << std::fixed << (minmax_coordinates.first)(1)*1e9 << " , " << (minmax_coordinates.second)(1)*1e9 << ") [nm]\n";
		std::cout << "    z (" << std::fixed << (minmax_coordinates.first)(2)*1e9 << " , " << (minmax_coordinates.second)(2)*1e9 << ") [nm]\n";

		std::cout.flags(f);

		return minmax_coordinates;
	};

	// step the simulation in time
	void step(double dt)
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
	void population_profiler(const unsigned& max_history, std::fstream& file, std::fstream& debug_file)
	{

		if (_population_probe.history < max_history)
		{
			// open the debug file
			if (! debug_file.is_open())
			{
				debug_file.open(_output_directory.path() / "debug.dat", std::ios::out);
				file << std::showpos << std::scientific;
			}

			t_region& bulk = _regions[1];

			_population_probe.history += 1;

			unsigned i;
			for (const auto& m_particle_ptr: bulk.particles())
			{
				i = std::floor((m_particle_ptr->pos(_population_probe.dim) - bulk.lower_corner(_population_probe.dim))/_population_probe.dL);
				if ((i >= _population_probe.number_of_sections) or (i<0))
				{
					debug_file << "index out of bound when making population profile\n";
					debug_file << "particle position = " << m_particle_ptr->pos(1) << " , bulk limits = [ " << bulk.lower_corner(_population_probe.dim) << " , " << bulk.upper_corner(1) << " ]" << std::endl;
				}
				else
				{
					_population_probe.profile[i] += 1;
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
				for (unsigned i=0; i < _population_probe.number_of_sections; ++i)
				{
					file << double(i)*_population_probe.dL << " ";
				}
				file << std::endl;
			}

			file << time() << " ";
			for (auto& element: _population_probe.profile)
			{
				file << double(element)/double(_population_probe.history)/_population_probe.dV << " ";
				element = 0;
			}
			file << std::endl;

			_population_probe.history = 0;
		}
	};

	// save the net currents for each region calculated by counting in and out flow of particles in each contact
	void save_region_current(const unsigned& max_history, std::fstream& current_file, const double& time_step)
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
			double cross_section = (bulk.upper_corner(0)-bulk.lower_corner(0)) * (bulk.upper_corner(2)-bulk.lower_corner(2));
			double elapsed_time = double(_history_of_region_currents)*time_step;

			current_file << time() << " ";
			for (auto&& m_region: _regions)
			{
				current_file << double(m_region.particle_flow())/elapsed_time/cross_section;
				current_file << " ";
				m_region.reset_particle_flow();
			}
			current_file << std::endl;

			_history_of_region_currents = 0;
		}
	};
	
	// read in the coordinate of all the cnt segments or molecules and create the scatterer objects that manage particle hopping between the sites
	void create_scatterers_without_orientation(const std::experimental::filesystem::path& input_path, const std::experimental::filesystem::path& output_path)
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
						arma::vec pos;
						std::string token;
						int i = 0;
						while(std::getline(iss,token,','))
						{
							pos(i) = 1.e-9*std::stod(token);
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

	// read in the coordinate of all the cnt segments or molecules and create the scatterer objects that manage particle hopping between the sites
	void create_scatterers_with_orientation(const std::experimental::filesystem::path& input_path, const std::experimental::filesystem::path& output_path)
	{

		std::cout << "\n\n...you are using the new_create_scatterers function!!!...\n\n";

		std::cout << "this is the input path: " << input_path << std::endl;
		std::ifstream file;
		std::string line;

		int num = 1;
		int max_num = 5;
		std::string base = "tube";
		std::string extension = ".dat";
		std::string filename = input_path / (base+std::to_string(num)+extension);

		file.open(filename);
		std::regex tube_rgx("tube");

		// convert a long string in the form of " text ; num0 , num1 , num2 ; num0 , num1 , num2 ;text" into a list of strings with form "num0 , num1 , num2"
		auto get_nodes = [](std::string str) -> std::list<std::string> {
			std::list<std::string> output;
			std::istringstream iss(str);
			std::string token;
			while (std::getline(iss, token, ';'))
			{
				int count = 0;
				for (auto t : token)
				{
					if (t == ',')	count ++;
				}
				if (count == 2)
				{
					output.emplace_back(token);
				}
			}
			return output;
		};

		// convert a string in the form of "num0 , num1 , num2" into an array of numbers
		auto get_position = [](std::string str) {
			int count = 0;
			for (auto s : str)
			{
				if (s == ',')	count ++;
			}
			arma::rowvec pos(3);
			if (count == 2)
			{
				std::istringstream iss(str);
				std::string token;
				std::getline(iss,token,',');
				pos(0) = std::stod(token);
				std::getline(iss,token,',');
				pos(1) = std::stod(token);
				std::getline(iss,token);
				pos(2) = std::stod(token);
			}
			pos = pos*1.e-9;
			return pos;
		};

		// loop over files
		while ((file.is_open()) and (num<max_num))
		{
			std::cout << "reading data from file: \"" << filename << "\"\r" << std::flush;

			while(std::getline(file, line, '\n'))
			{

				// find the new tube coordinates by finding the keywork "tube"
				if (std::regex_search(line,tube_rgx))
				{
					auto nodes = get_nodes(line);

					arma::mat tube_coordinates(nodes.size(),3);
					unsigned i=0;
					for (const auto& node: nodes)
					{
						tube_coordinates.row(i) = get_position(node);
						i++;
					}

					for (unsigned i=0; i<tube_coordinates.n_rows; i++)
					{
						arma::rowvec pos1;
						arma::rowvec pos2;
						if (i==0)
						{
							pos1 = tube_coordinates.row(i);
							pos2 = tube_coordinates.row(i+1);
						}
						else if(i==(tube_coordinates.n_rows-1))
						{
							pos1 = tube_coordinates.row(i-1);
							pos2 = tube_coordinates.row(i);
						}
						else
						{
							pos1 = tube_coordinates.row(i-1);
							pos2 = tube_coordinates.row(i+1);
						}

						arma::rowvec orientation = arma::normalise(pos2-pos1);

						_all_scat_list.push_back(std::make_shared<mc::discrete_forster_scatter>());
						_all_scat_list.back()->set_pos(tube_coordinates.row(i));
						_all_scat_list.back()->set_orientation(orientation);
					}

				}
			}

			file.close();
			num++;
			filename = input_path / (base+std::to_string(num)+extension);
			file.open(filename);

		}

		std::cout << "\n\ntotal number of discrete_forster_scatterer: " << _all_scat_list.size() << std::endl;
	};

  // high level method to calculate proper scattering table
	void initialize_scattering_table();

	// method to calculate scattering rate via forster method
	scattering_struct create_forster_scatt_table(double gamma_0, double r_0);

	// method to calculate scattering rate via davoody et al. method
	scattering_struct create_davoody_scatt_table(const cnt& d_cnt, const cnt& a_cnt);


}; // end class discrete_forster_monte_carlo

} // end namespace mc

#endif // discrete_forster_monte_carlo_h
