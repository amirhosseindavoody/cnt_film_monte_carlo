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

#include "../../lib/json.hpp"

#include "../exciton_transfer/cnt.h"
#include "../exciton_transfer/exciton_transfer.h"
#include "./region.h"
#include "./particle.h"
#include "./scatterer.h"
#include "./scattering_struct.h"

namespace mc
{

class discrete_forster_monte_carlo
{

private:
  struct population_profile {
    std::vector<unsigned> profile;  // number of particles in each segment
    unsigned number_of_sections=0;    // number of sections in the profiler
    unsigned history=0;               // history of the population_profiler
    double dL;     // length of each segment in population profiler
    double dV;     // volume of each segment in population profiler
    unsigned dim;  // this is the dimension along which the population profiler
                    // is working
  };

  double _time;
  double _max_hopping_radius;

  nlohmann::json _json_prop;  // input properties of the whole mc simulation in json format

  std::array<region_class, 3> _regions;  // simulation regions: contacts and bulk. Since
                                         // discrete_forster_region data type has unique_ptr as data
                                         // members, we can't use vectors and dynamically push_back new
                                         // regions. If we want to do that we need to implement the move
                                         // and copy constructors explicitly.

  std::vector<scatterer> _all_scat_list;


  std::pair<arma::vec, arma::vec> _domain;

  // number of particles in the contacts
  unsigned _number_of_contact1_particles, _number_of_contact2_particles;

  // this is the address of the output_directory and input_directory
  std::experimental::filesystem::directory_entry _output_directory, _input_directory;

  unsigned _history_of_region_currents = 0;  // this is the number of steps that the net
                                             // _current in the regions have‍‍‍‍ been recorded
  population_profile _population_probe;  // this is the population profile of particles through
                                         // the simulation domain along the z-axis

  scattering_struct _scat_table;  // instantiation of the scattering table for
                                  // discrete mesh points

  std::vector<std::vector<scatterer*>> _scat_buckets;  // pointers to scatterers to divide the scatterers
                                                       // into multiple buckets based on their position in
                                                       // space

  std::vector<particle> _particle_list;

  free_flight _ff;

public:
	// default constructor
	discrete_forster_monte_carlo()=delete;
	
	// constructure with json input file
	discrete_forster_monte_carlo(const nlohmann::json& j) {
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

		// set maximum hopping radius
    _max_hopping_radius = double(j["max hopping radius [m]"]);
		std::cout << "maximum hopping radius: " << _max_hopping_radius*1.e9 << " [nm]\n";

	};

	// get the mc simulation time
	const double& time() const {
		return _time;
	};

	// increase simulation time by dt
	void increase_time(const double& dt) {
		_time += dt;
	};

	// get constant reference to the output_directory
	const std::experimental::filesystem::directory_entry& output_directory() const{
		return _output_directory;
	};

	// get constant reference to the input_directory
	const std::experimental::filesystem::directory_entry& input_directory() const{
		return _input_directory;
	};

	// get constant reference to the output_directory
	const std::experimental::filesystem::path& output_path() const{
		return _output_directory.path();
	};

	// get constant reference to the input_directory
	const std::experimental::filesystem::path& input_path() const {
		return _input_directory.path();
	};

	// returns the number of particles
	unsigned number_of_particles() const {
		unsigned number_of_particles = 0;
		for (const auto& m_region: _regions)
		{
			number_of_particles += m_region.number_of_particles();
		}
		return number_of_particles;
	};

	//initialize the simulation condition
	void init() {

    _scat_table = create_scattering_table(_json_prop);
    _all_scat_list = create_scatterers(_json_prop);

    set_scat_table(_scat_table, _all_scat_list);

    _domain = find_simulation_domain();

    create_scatterer_buckets(_domain, _max_hopping_radius, _all_scat_list, _scat_buckets);
    // set_max_rate(_max_hopping_radius, _all_scat_list);

    _number_of_contact1_particles = 1100;
    _number_of_contact2_particles = 0;

    std::cout << "ready to create particles...\n";
    std::cin.ignore();

    _ff = free_flight();

    _particle_list =
        create_particles(_domain, _all_scat_list, _number_of_contact1_particles, _number_of_contact2_particles, &_ff);

    std::exit(0);

    double x_min = (_domain.first)(0);
    double x_max = (_domain.second)(0);

    double y_min = (_domain.first)(1);
    double y_1 = (_domain.first)(1) + 0.1 * ((_domain.second)(1) - (_domain.first)(1));
    double y_2 = (_domain.first)(1) + 0.9 * ((_domain.second)(1) - (_domain.first)(1));
    double y_max = (_domain.second)(1);

    double z_min = (_domain.first)(2);
    double z_max = (_domain.second)(2);

    arma::vec contact_1_lower_corner = {x_min, y_min, z_min};
    arma::vec contact_1_upper_corner = {x_max, y_1, z_max};
    arma::vec bulk_lower_corner = {x_min, y_1, z_min};
    arma::vec bulk_upper_corner = {x_max, y_2, z_max};
    arma::vec contact_2_lower_corner = {x_min, y_2, z_min};
    arma::vec contact_2_upper_corner = {x_max, y_max, z_max};

    _regions[0].set_borders(contact_1_lower_corner, contact_1_upper_corner);
    _regions[1].set_borders(bulk_lower_corner, bulk_upper_corner);
    _regions[2].set_borders(contact_2_lower_corner, contact_2_upper_corner);

    for (auto& m_region : _regions) {
      m_region.create_scatterer_vector(_all_scat_list);
      m_region.create_pilot_list();
		}


		_number_of_contact1_particles = 1100;
		_number_of_contact2_particles = 0;

		std::cout << "\n"
              << "number of particles in the contact 1: " << _number_of_contact1_particles << "\n"
							<< "number of particles in the contact 2: " << _number_of_contact2_particles << "\n"
              << "\n";

		_regions.front().populate(_number_of_contact1_particles);
		_regions.back().populate(_number_of_contact2_particles);

		region_class* bulk = &(_regions[1]);

		_population_probe.number_of_sections = 10;
		_population_probe.dim = 1;
		_population_probe.profile = std::vector<unsigned>(_population_probe.number_of_sections, 0);
		_population_probe.history = 0;
		_population_probe.dL = (bulk->upper_corner(_population_probe.dim)-bulk->lower_corner(_population_probe.dim))/double(_population_probe.number_of_sections);
		_population_probe.dV = bulk->volume()/(_population_probe.number_of_sections);

	};

  // high level function to create scatterers vector based on the json inpu
  std::vector<scatterer> create_scatterers(nlohmann::json j){
    if (j.count("mesh type") == 0) throw std::invalid_argument("specify mesh type in input json");

    std::vector<scatterer> scat_list;

    std::string mesh_type = _json_prop["mesh type"];
    if (_json_prop["mesh type"] == "crystalline") {
      scat_list = create_crystalline_structure();
    } else if (_json_prop["mesh type"] == "realistic") {
      create_scatterers_with_orientation(input_directory().path(), output_directory().path());
    } else if (_json_prop["mesh type"] == "fiber") {
      scat_list = create_scatterer_from_fiber(input_directory().path());
    } else {
      throw std::invalid_argument("invalid mesh type in input json");
    }

    std::cout << "\n\n"
              << "total number of scatterers: " << scat_list.size() << std::endl;
    return scat_list;
  }

  // create particles with a linear density profile in y direction
  std::vector<particle> create_particles( const std::pair<arma::vec, arma::vec>& domain,
      const std::vector<scatterer>& scat_list, int left_pop, int right_pop, const free_flight* ff) {

    std::cout << "\n"
              << "create particles list: ";

    std::vector<particle> p_list;
    int no_of_sections = 10;

    double y_min = domain.first(1);
    double y_max = domain.second(1);

    for (int i=0; i<no_of_sections; ++i){

      std::cout << i << ",";

      double dp = double(right_pop - left_pop) / double(no_of_sections - 1);
      int n_particle = std::round(left_pop + i * dp);

      double dy = (y_max-y_min)/no_of_sections;
      double y1 = y_min+i*dy;
      double y2 = y_min+(i+1)*dy;

      std::vector<const scatterer*> s_list;
      for (const scatterer& s: scat_list){
        if (y1<=s.pos(1) && s.pos(1)<y2){
          s_list.emplace_back(&s);
        }
      }

      for (int n=0; n<n_particle; n++){
        int dice = std::rand()%s_list.size();
        const scatterer* s = s_list[dice];
        arma::vec pos = s->pos();
        p_list.push_back(particle(pos,ff,s));
      }
    }

    std::cout << "...done!!!" << std::endl;

    return p_list;
  }

	// save the json properties that is read and parsed from the input_json file.
	void save_json_properties() {
		std::ofstream json_file;
		json_file.open(_output_directory.path() / "input.json", std::ios::out);
		json_file << std::setw(4) << _json_prop << std::endl;
		json_file.close();
	};

  // find minimum of the minimum coordinates of the scattering objects,
  // this function will effectively give us the simulation domain
	std::pair<arma::vec, arma::vec> find_simulation_domain() {
		
    arma::vec min_coor = _all_scat_list.front().pos();
    arma::vec max_coor = _all_scat_list.front().pos();

    std::pair<arma::vec, arma::vec> minmax_coordinates(_all_scat_list.front().pos(),
																											 _all_scat_list.front().pos());

    for (const auto& s : _all_scat_list) {
      for (int i = 0; i < 3; ++i) {
        min_coor(i) = min_coor(i) > s.pos(i) ? s.pos(i) : min_coor(i);
        max_coor(i) = max_coor(i) < s.pos(i) ? s.pos(i) : max_coor(i);
      }
    }

    std::ios::fmtflags f(
        std::cout.flags());  // save cout flags to be reset after printing

    std::cout << "\n simulation domain:\n";
    std::cout << "    x (" << std::fixed << std::showpos << min_coor(0) * 1e9
              << " , " << max_coor(0) * 1e9 << ") [nm]\n";
    std::cout << "    y (" << std::fixed << min_coor(1) * 1e9 << " , "
              << max_coor(1) * 1e9 << ") [nm]\n";
    std::cout << "    z (" << std::fixed << min_coor(2) * 1e9 << " , "
              << max_coor(2) * 1e9 << ") [nm]\n";

    std::cout.flags(f);  // reset the cout flags

    return {min_coor, max_coor};
	};

	// step the simulation in time
	void step(double dt) {

		for (auto&& p=_regions[0].particles().begin(); p!=_regions[0].particles().end(); ++p) {
			(*p)->step(dt,_domain, _max_hopping_radius);
			_regions[1].enlist(p,&(_regions[0]));
		}


		for (auto&& p=_regions[2].particles().begin(); p!=_regions[2].particles().end(); ++p) {
      (*p)->step(dt, _domain, _max_hopping_radius);
      _regions[1].enlist(p, &(_regions[2]));
		}

		for (auto&& p=_regions[1].particles().begin(); p!=_regions[1].particles().end(); ++p) {
			(*p)->step(dt, _domain, _max_hopping_radius);
			if (not _regions[2].enlist(p,&(_regions[1]))) {
				_regions[0].enlist(p,&(_regions[1]));
			}
		}

		// dump the new particles into the particles list in each region
		for (auto& m_region: _regions) {
			m_region.dump_new_particles();
		}

		increase_time(dt);

	};

	// repopulate contacts
	void repopulate_contacts() {
		_regions[0].populate(_number_of_contact1_particles);
		_regions[1].populate(_number_of_contact2_particles);
	};

	// calculate and save the population profile
	void population_profiler(const unsigned& max_history, std::fstream& file, std::fstream& debug_file) {

    // open the debug file
    if (! debug_file.is_open()) {
      debug_file.open(_output_directory.path() / "debug.dat", std::ios::out);
      file << std::showpos << std::scientific;
    }

    region_class* bulk = &(_regions[1]);

    _population_probe.history++;

    unsigned i;
    for (const auto& p: bulk->particles()) {
      i = std::floor((p->pos(_population_probe.dim) - bulk->lower_corner(_population_probe.dim))/_population_probe.dL);
      if ((i >= _population_probe.number_of_sections) or (i<0)) {
        debug_file << "index out of bound when making population profile\n";
        debug_file << "particle position = " << p->pos(_population_probe.dim) << " , bulk limits = [ " << bulk->lower_corner(_population_probe.dim) << " , " << bulk->upper_corner(_population_probe.dim) << " ]" << std::endl;

        std::cout << "_population_probe.number_of_sections = " << _population_probe.number_of_sections << "\ni = " << i << "\n";

        std::cout << "particle position = " << p->pos(_population_probe.dim) << "\n"
                  << "bulk limits = [" << bulk->lower_corner(_population_probe.dim) << " , " << bulk->upper_corner(_population_probe.dim) << "]\n"
                  << "bulk limits - particle position = [" << bulk->lower_corner(_population_probe.dim) - p->pos(_population_probe.dim) << " , " << bulk->upper_corner(_population_probe.dim) - p->pos(_population_probe.dim) << "]\n";
        throw std::range_error("particle is out of bound when calculating population profile");

      } else {
        _population_probe.profile[i]++;
      }
    }

    if (_population_probe.history == max_history) {
      if (!file.is_open()) {
        file.open(_output_directory.path() / "population_profile.dat",
                  std::ios::out);
        file << std::showpos << std::scientific;

        // store the position of the middle point of each section
        file << time() << " ";
        for (unsigned i = 0; i < _population_probe.number_of_sections; ++i) {
          file << double(i) * _population_probe.dL << " ";
        }
        file << std::endl;
      }

      file << time() << " ";
      for (auto& element : _population_probe.profile) {
        file << double(element) / double(_population_probe.history) / _population_probe.dV << " ";
        element = 0;
      }
      file << std::endl;

      _population_probe.history = 0;
    }
	};

	// save the net currents for each region calculated by counting in and out flow of particles in each contact
	void save_region_current(const unsigned& max_history, std::fstream& current_file, const double& time_step) {

    _history_of_region_currents++;

		if (_history_of_region_currents%max_history==0) {
			if (! current_file.is_open()) {
				current_file.open(_output_directory.path() / "region_current.dat", std::ios::out);
				current_file << std::showpos << std::scientific;
			}

			region_class* bulk = &(_regions[1]);
			double cross_section = (bulk->upper_corner(0)-bulk->lower_corner(0)) * (bulk->upper_corner(2)-bulk->lower_corner(2));
			double elapsed_time = double(_history_of_region_currents)*time_step;

			current_file << time() << " ";
			for (auto&& m_region: _regions) {
				current_file << double(m_region.particle_flow())/elapsed_time/cross_section;
				current_file << " ";
				m_region.reset_particle_flow();
			}
			current_file << std::endl;
		}
	};

  // read in the coordinate of all the cnt segments or molecules and
  // create the scatterer objects that manage particle hopping between the
  // sites
  void create_scatterers_with_orientation(
      const std::experimental::filesystem::path& input_path,
      const std::experimental::filesystem::path& output_path) {
    std::cout << "this is the input path: " << input_path << std::endl;
    std::ifstream file;
    std::string line;

    int num = 1;
    int max_num = 5;  // maximum number of files to read
    std::string base = "tube";
    std::string extension = ".dat";
    std::string filename =
        input_path / (base + std::to_string(num) + extension);

    file.open(filename);
    std::regex tube_rgx("tube");

    // convert a long string in the form of " text ; num0 , num1 , num2 ;
    // num0 , num1 , num2 ;text" into a list of strings with form "num0 ,
    // num1 , num2"
    auto get_nodes = [](std::string str) -> std::list<std::string> {
      std::list<std::string> output;
      std::istringstream iss(str);
      std::string token;
      while (std::getline(iss, token, ';')) {
        int count = 0;
        for (auto t : token) {
          if (t == ',') count++;
        }
        if (count == 2) {
          output.emplace_back(token);
        }
      }
      return output;
    };

    // convert a string in the form of "num0 , num1 , num2" into an array
    // of numbers
    auto get_position = [](std::string str) {
      int count = 0;
      for (auto s : str) {
        if (s == ',') count++;
      }
      arma::rowvec pos(3);
      if (count == 2) {
        std::istringstream iss(str);
        std::string token;
        std::getline(iss, token, ',');
        pos(0) = std::stod(token);
        std::getline(iss, token, ',');
        pos(1) = std::stod(token);
        std::getline(iss, token);
        pos(2) = std::stod(token);
      }
      pos = pos * 1.e-9;
      return pos;
    };

    // loop over files
    while ((file.is_open()) and (num <= max_num)) {
      std::cout << "reading data from file: \"" << filename << "\"\r"
                << std::flush;

      while (std::getline(file, line, '\n')) {
        // find the new tube coordinates by finding the keywork "tube"
        if (std::regex_search(line, tube_rgx)) {
          auto nodes = get_nodes(line);

          arma::mat tube_coordinates(nodes.size(), 3);
          unsigned i = 0;
          for (const auto& node : nodes) {
            tube_coordinates.row(i) = get_position(node);
            i++;
          }

          for (unsigned i = 0; i < tube_coordinates.n_rows; i++) {
            arma::rowvec pos1;
            arma::rowvec pos2;
            if (i == 0) {
              pos1 = tube_coordinates.row(i);
              pos2 = tube_coordinates.row(i + 1);
            } else if (i == (tube_coordinates.n_rows - 1)) {
              pos1 = tube_coordinates.row(i - 1);
              pos2 = tube_coordinates.row(i);
            } else {
              pos1 = tube_coordinates.row(i - 1);
              pos2 = tube_coordinates.row(i + 1);
            }
            arma::rowvec orientation = arma::normalise(pos2 - pos1);

            _all_scat_list.push_back(scatterer());
            _all_scat_list.back().set_pos(tube_coordinates.row(i).t());
            _all_scat_list.back().set_orientation(orientation.t());
          }
        }
      }

      file.close();
      num++;
      filename = input_path / (base + std::to_string(num) + extension);
      file.open(filename);
    }
  };

  // create a crystalline mesh structure
  std::vector<scatterer> create_crystalline_structure();

  // high level method to calculate proper scattering table
	scattering_struct create_scattering_table(nlohmann::json j);

	// method to calculate scattering rate via forster method
	scattering_struct create_forster_scatt_table(double gamma_0, double r_0);

	// method to calculate scattering rate via davoody et al. method
	scattering_struct create_davoody_scatt_table(const cnt& d_cnt, const cnt& a_cnt);

  // read in the coordinate of all the cnt segments or molecules and create the scatterer objects that manage
  // particle hopping between the sites
  std::vector<scatterer> create_scatterer_from_fiber(const std::experimental::filesystem::path& input_path) {
    std::cout << "\n"
              << "create scatterers in fiber sctructure:\n"
              << "input path " << input_path << std::endl;

    std::ifstream pos_file;
    pos_file.open(input_path / "single_cnt.pos.dat");
    std::ifstream orient_file;
    orient_file.open(input_path / "single_cnt.orient.dat");

    arma::mat pos;
    pos.load(pos_file);
    pos *= 1.e-9;

    arma::mat orient;
    orient.load(orient_file);

    pos_file.close();
    orient_file.close();

    std::vector<scatterer> scat_list;

    scat_list.resize(pos.n_rows);

    for (unsigned i = 0; i < pos.n_rows; ++i) {
      scat_list[i].set_pos(pos.row(i).t());
      scat_list[i].set_orientation(orient.row(i).t());
    }

    std::cout << "done!!!" << std::endl;

    return scat_list;
  };

  // divide scatterers into buckets based on their location, and set the
  // pointers to enclosing and neighboring buckets for each scatterer object
  void create_scatterer_buckets(const std::pair<arma::vec, arma::vec> domain, const double max_hopping_radius, std::vector<scatterer>& scat_list,
                                std::vector<std::vector<scatterer*>>& scat_buckets) {
    using namespace std;
    
    std::cout << "\n" 
              << "finding scatterer buckets: ";

    double xmin = (domain.first)(0);
    double xmax = (domain.second)(0);
    int nx = std::ceil((xmax - xmin) / max_hopping_radius) + 1;

    double ymin = (domain.first)(1);
    double ymax = (domain.second)(1);
    int ny = std::ceil((ymax - ymin) / max_hopping_radius) + 1;

    double zmin = (domain.first)(2);
    double zmax = (domain.second)(2);
    int nz = std::ceil((zmax - zmin) / max_hopping_radius) + 1;

    scat_buckets.resize(nx*ny*nz);

    for (scatterer& s : scat_list) {
      int ix = (s.pos(0) - xmin) / max_hopping_radius;
      int iy = (s.pos(1) - ymin) / max_hopping_radius;
      int iz = (s.pos(2) - zmin) / max_hopping_radius;
      int idx = ix + iy * nx + iz * nx * ny;
      scat_buckets[idx].push_back(&s);
    }

    for (scatterer& s : scat_list) {
      int ix = (s.pos(0) - xmin) / max_hopping_radius;
      int iy = (s.pos(1) - ymin) / max_hopping_radius;
      int iz = (s.pos(2) - zmin) / max_hopping_radius;

      for (int i : {ix - 1, ix, ix + 1}) {
        for (int j : {iy - 1, iy, iy + 1}) {
          for (int k : {iz - 1, iz, iz + 1}) {
            if (i > -1 && i < nx && j > -1 && j < ny && k > -1 && k < nz) {
              unsigned idx = i + j * nx + k * nx * ny;
              s.close_scats.push_back(&(scat_buckets[idx]));
            }
          }
        }
      }
    }

    std::cout << "done!\n";

  }

  // set the pointer to scattering table struct for all scatterer objects
  void set_scat_table(const scattering_struct& scat_tab,
                      std::vector<scatterer>& scat_list) {
    for (auto& s : scat_list) {
      s.scat_tab = &scat_tab;
    }
  }

  // set the max scattering rate for all the scatterers
  void set_max_rate(const double max_hopping_radius, std::vector<scatterer>& scat_list){
    std::cout << "\nsetting max rate in scatterers:" << std::endl;
    progress_bar prog(scat_list.size(), "setting max rate in scatterers");
    
    for (auto& s: scat_list){
      s.set_max_rate(max_hopping_radius);
      prog.step();
    }
  }

}; // end class discrete_forster_monte_carlo

} // end namespace mc

#endif // discrete_forster_monte_carlo_h
