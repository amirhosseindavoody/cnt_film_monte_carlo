#ifndef monte_carlo_h
#define monte_carlo_h

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

class monte_carlo
{

private:
  typedef std::experimental::filesystem::path            path_t;
  typedef std::experimental::filesystem::directory_entry directory_t;
  typedef std::pair<arma::vec, arma::vec>                domain_t;
  typedef std::vector<std::vector<scatterer*>>           bucket_t;

  // elapsed simulation time
  double _time;

  // maximum hopping radius considered in the simulation
  double _max_hopping_radius;

  // input properties of the whole mc simulation in json format
  nlohmann::json _json_prop;  

  // list of all scatterer object in the simulation
  std::vector<scatterer> _all_scat_list;

  // minimum and maximum coordinates of the simulation domain
  domain_t _domain;

  // number of particles in the contacts
  unsigned _c1_pop, _c2_pop;

  // this is the address of the output_directory and input_directory
  directory_t _output_directory, _input_directory;

  // instantiation of the scattering table for discrete mesh points
  scattering_struct _scat_table;

  // pointers to scatterers to divide the scatterers into multiple buckets based on their position in space
  bucket_t _scat_buckets;  

  // pointers to scatterers in contact 1 and 2
  std::vector<const scatterer*> _c1_scat, _c2_scat;

  // number of segments that defines contacts
  unsigned _n_seg=0;

  // list of all particles in the simulation
  std::vector<particle> _particle_list;

  // pilot object for the free flight section of the monte carlo
  free_flight _ff;

  // file objects for saving population profile and current data
  std::fstream _pop_file, _curr_file;

 public:
  // default constructor
  monte_carlo() = delete;

  // constructure with json input file
  monte_carlo(const nlohmann::json& j) {

    std::cout << "\n"
              << "ready properties from json file"
              << "\n";

    // store the json properties for use in other methods
    _json_prop = j;

    // set the output directory
    std::string directory_path = j["output directory"];
    bool        keep_old_data = true;
    if (j.count("keep old results") == 1) {
      keep_old_data = j["keep old results"];
    }
    _output_directory = prepare_directory(directory_path, keep_old_data);

    // set the input directory for mesh information
    directory_path = j["mesh input directory"];
    _input_directory = check_directory(directory_path, false);

    // set maximum hopping radius
    _max_hopping_radius = double(j["max hopping radius [m]"]);
		std::cout << "maximum hopping radius: " << _max_hopping_radius*1.e9 << " [nm]\n";

    // set the number of segments along y axis
    _n_seg = j["number of segments"];
    std::cout << "number of segments: " << _n_seg << std::endl;

	};

	// get the mc simulation time
  const double& time() const { return _time; };

  // get constant reference to the output_directory
  const directory_t& output_directory() const { return _output_directory; };

  // get constant reference to the input_directory
  const directory_t& input_directory() const { return _input_directory; };

  // get constant reference to the output_directory
  const path_t& output_path() const { return _output_directory.path(); };

  // get constant reference to the input_directory
  const path_t& input_path() const { return _input_directory.path(); };

  // returns the number of particles
  unsigned number_of_particles() const { return _particle_list.size(); };

  // initialize the simulation condition
  void init() {

    _scat_table = create_scattering_table(_json_prop);
    _all_scat_list = create_scatterers(_json_prop);

    set_scat_table(_scat_table, _all_scat_list);

    _domain = find_simulation_domain();

    create_scatterer_buckets(_domain, _max_hopping_radius, _all_scat_list, _scat_buckets);
    set_max_rate(_max_hopping_radius, _all_scat_list);

    _c1_scat = contact_scats(_all_scat_list, _n_seg, 1, _domain);
    _c2_scat = contact_scats(_all_scat_list, _n_seg, _n_seg, _domain);

    _c1_pop = 1100;
    _c2_pop = 0;

    _ff = free_flight();

    _particle_list = create_particles(_domain, _n_seg, _all_scat_list, _c1_pop, _c2_pop, &_ff);

    _pop_file.open(_output_directory.path() / "population_profile.dat", std::ios::out);
    _curr_file.open(_output_directory.path() / "region_current.dat", std::ios::out);
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
  std::vector<particle> create_particles( const domain_t& domain, const unsigned n_seg,
      const std::vector<scatterer>& scat_list, int left_pop, int right_pop, const free_flight* ff) {

    std::cout << "\n"
              << "create particles list:...";

    std::vector<particle> p_list;

    double y_min = domain.first(1);
    double y_max = domain.second(1);

    double dy = (y_max - y_min) / double(n_seg);
    double dp = double(right_pop - left_pop) / (double(n_seg) - 1);

    for (unsigned i=0; i<n_seg; ++i){


      int n_particle = std::round(left_pop + double(i) * dp);

      std::cout << "("<< i << "," << n_particle << ") ,";
      
      double y1 = y_min + double(i) * dy;
      double y2 = y1 + dy;

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

  // find minimum of the minimum coordinates of the scattering objects, this function will effectively give us the
  // simulation domain
  domain_t find_simulation_domain() {
    arma::vec min_coor = _all_scat_list.front().pos();
    arma::vec max_coor = _all_scat_list.front().pos();

    for (const auto& s : _all_scat_list) {
      for (int i = 0; i < 3; ++i) {
        min_coor(i) = min_coor(i) > s.pos(i) ? s.pos(i) : min_coor(i);
        max_coor(i) = max_coor(i) < s.pos(i) ? s.pos(i) : max_coor(i);
      }
    }

    std::ios::fmtflags f(std::cout.flags());  // save cout flags to be reset after printing

    std::cout << std::fixed << std::showpos;

    std::cout << "\n simulation domain:\n";
    std::cout << "    x (" << min_coor(0) * 1e9 << " , " << max_coor(0) * 1e9 << ") [nm]\n";
    std::cout << "    y (" << min_coor(1) * 1e9 << " , " << max_coor(1) * 1e9 << ") [nm]\n";
    std::cout << "    z (" << min_coor(2) * 1e9 << " , " << max_coor(2) * 1e9 << ") [nm]\n";

    std::cout.flags(f);  // reset the cout flags

    return {min_coor, max_coor};
  };

  // step the simulation in time
	void step(double dt) {

    for (auto& p: _particle_list){
      p.step(dt, _domain, _max_hopping_radius);
    }

    // increase simulation time
    _time += dt;
	};

  // read in the coordinate of all the cnt segments or molecules and create the scatterer objects that manage
  // particle hopping between the sites
  
  void create_scatterers_with_orientation(const path_t& input_path, const path_t& output_path) {
    std::cout << "this is the input path: " << input_path << std::endl;
    std::ifstream file;
    std::string   line;

    int         num = 1;
    int         max_num = 5;  // maximum number of files to read
    std::string base = "tube";
    std::string extension = ".dat";
    std::string filename = input_path / (base + std::to_string(num) + extension);

    file.open(filename);
    std::regex tube_rgx("tube");

    // convert a long string in the form of " text ; num0 , num1 , num2 ;
    // num0 , num1 , num2 ;text" into a list of strings with form "num0 ,
    // num1 , num2"
    auto get_nodes = [](std::string str) -> std::list<std::string> {
      std::list<std::string> output;
      std::istringstream     iss(str);
      std::string            token;
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
        std::string        token;
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
      std::cout << "reading data from file: \"" << filename << "\"\r" << std::flush;

      while (std::getline(file, line, '\n')) {
        // find the new tube coordinates by finding the keywork "tube"
        if (std::regex_search(line, tube_rgx)) {
          auto nodes = get_nodes(line);

          arma::mat tube_coordinates(nodes.size(), 3);
          unsigned  i = 0;
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
  std::vector<scatterer> create_scatterer_from_fiber(const path_t& input_path) {
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

    for (unsigned i = 0; i < scat_list.size(); ++i) {
      scat_list[i].set_pos(pos.row(i).t());
      scat_list[i].set_orientation(orient.row(i).t());
    }

    std::cout << "done!!!" << std::endl;

    return scat_list;
  };

  // divide scatterers into buckets based on their location, and set the pointers to enclosing and neighboring buckets
  // for each scatterer object
  void create_scatterer_buckets(const domain_t domain, const double radius, std::vector<scatterer>& scat_list,
                                bucket_t& scat_buckets) {
    using namespace std;
    
    std::cout << "\n" 
              << "finding scatterer buckets: ";

    double xmin = (domain.first)(0);
    double xmax = (domain.second)(0);
    int nx = std::ceil((xmax - xmin) / radius) + 1;

    double ymin = (domain.first)(1);
    double ymax = (domain.second)(1);
    int    ny = std::ceil((ymax - ymin) / radius) + 1;

    double zmin = (domain.first)(2);
    double zmax = (domain.second)(2);
    int    nz = std::ceil((zmax - zmin) / radius) + 1;

    scat_buckets.resize(nx*ny*nz);

    for (scatterer& s : scat_list) {
      int ix = (s.pos(0) - xmin) / radius;
      int iy = (s.pos(1) - ymin) / radius;
      int iz = (s.pos(2) - zmin) / radius;
      int idx = ix + iy * nx + iz * nx * ny;
      scat_buckets[idx].push_back(&s);
    }

    for (scatterer& s : scat_list) {
      int ix = (s.pos(0) - xmin) / radius;
      int iy = (s.pos(1) - ymin) / radius;
      int iz = (s.pos(2) - zmin) / radius;

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
  void set_scat_table(const scattering_struct& scat_tab, std::vector<scatterer>& scat_list) {
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

  // repopulate contacts
  void repopulate_contacts() {
    double ymin = _domain.first(1);
    double ymax = _domain.second(1);
    double dy = (ymax - ymin) / double(_n_seg);

    double y1 = ymin;
    double y2 = ymin + dy;
    repopulate(y1, y2, _c1_pop, &_ff, _c1_scat, _particle_list);

    y1 = ymin + double(_n_seg - 1) * dy;
    y2 = ymax;
    repopulate(y1, y2, _c2_pop, &_ff, _c2_scat, _particle_list);
  };

  // take all the particles between ymin and ymax region and recycle them and populate the region with new particles
  void repopulate(const double ymin, const double ymax, const unsigned n_particle, const free_flight* ff,
                  const std::vector<const scatterer*>& s_list,
                  std::vector<particle>& p_list) {
    
    unsigned j=p_list.size();

    for (unsigned i = 0; i < j;) {
      if (p_list[i].pos(1) >= ymin && p_list[i].pos(1) <= ymax) {
        --j;
        std::swap(p_list[i], p_list[j]);
      } else {
        ++i;
      }
    }

    int dice=0;
    unsigned n=0;
    unsigned final_size = j+n_particle;

    unsigned j_lim = std::min(int(p_list.size()), int(final_size));
    
    for (;j < j_lim; ++j) {
      dice = std::rand() % s_list.size();
      p_list[j] = particle(s_list[dice]->pos(), ff, s_list[dice]);
      ++n;
    }

    for (;n<n_particle; ++n){
      dice = std::rand() % s_list.size();
      p_list.emplace_back(particle(s_list[dice]->pos(), ff, s_list[dice]));
    }

    p_list.resize(final_size);
  }

  // create a list of scatterer pointers in the contact number i
  std::vector<const scatterer*> contact_scats(const std::vector<scatterer>& s_list, const unsigned n_seg, const unsigned i,
                                        const domain_t& domain) {
    
    if (i==0 || i>n_seg)
      throw std::logic_error("input \"i\" must be between 1 and number of segments \"n_seg\"");
    
    double ymin = domain.first(1);
    double ymax = domain.second(1);
    double dy = (ymax-ymin)/double(n_seg);

    double y1 = ymin + double(i - 1) * dy;
    double y2 = ymin + double(i) * dy;

    std::vector<const scatterer*> c_list;

    for (auto& s: s_list){
      if (s.pos(1)>=y1 && s.pos(1)<=y2){
        c_list.push_back(&s);
      }
    }

    return c_list;
  }

  // calculate all the metrics needed from the experiment
  void save_metrics() {
    save_population_profile(_n_seg);
    save_currents(_n_seg);
  }

  // calculate and save population profile
  void save_population_profile(unsigned n) {
    std::vector<int> pop(n, 0);

    double ymax = (_domain.second)(1);
    double ymin = (_domain.first)(1);

    double dy = (ymax - ymin) / double(n);

    int i = 0;

    for (auto p : _particle_list) {
      i = (p.pos(1) - ymin) / dy;
      if (i > -1 && i < int(n)) {
        pop[i]++;
      }
    }

    if (!_pop_file.is_open()) {
      throw std::logic_error("Population file is not open: _pop_file !!!");
    }

    // std::cout << std::endl;
    // std::cout << "population profile:" << _time << "...";
    // for (unsigned i=0; i<pop.size(); ++i){
    //   std::cout << "(" << i << ","<< pop[i] << ") ,";
    // }
    // std::cout << std::endl;

    _pop_file << std::showpos << std::scientific << _time << " ";
    for (int& p : pop) {
      _pop_file << std::showpos << std::scientific << p << " ";
    }
    _pop_file << std::endl;
  }

  // calculate and save population profile
  void save_currents(unsigned n) {
    double ymax = (_domain.second)(1);
    double ymin = (_domain.first)(1);
    double dy = (ymax - ymin) / double(n);

    if (n < 4) {
      throw std::invalid_argument(
          "number of segments should be larger than 4 for the \"save_currents\" function: \"n\"");
    }

    std::vector<double> y = {ymin + 2 * dy, ymin + double(n - 2) * dy};

    std::vector<int> curr(2, 0);
    int              i = 0;

    for (i = 0; i < 2; ++i) {
      for (auto p : _particle_list) {
        if (p.old_pos(1) < y[i] && p.pos(1) >= y[i]) {
          curr[i]++;
        } else if (p.old_pos(1) < y[i] && p.pos(1) >= y[i]) {
          curr[i]--;
        }
      }
    }

    if (!_curr_file.is_open()) {
      throw std::logic_error("Current file is not open: _curr_file !!!");
    }

    _curr_file << std::showpos << std::scientific << _time << " ";
    for (auto& c : curr) {
      _curr_file << std::showpos << std::scientific << c << " ";
    }
    _curr_file << std::endl;
  }

  void get_scatterer_distribution(){
    double ymin = _domain.first(1);
    double ymax = _domain.second(1);
    double dy = (ymax - ymin) / double(_n_seg);

    std::vector<long> pop(_n_seg, 0);

    for (auto& s:_all_scat_list){
      int i = int(std::abs(s.pos(1) - ymin) / dy)%_n_seg;
      pop[i]++;
    }

    std::vector<double> pos(_n_seg, 0);
    for (unsigned i=0; i<pos.size(); ++i){
      pos[i] = ymin + (double(i) + 0.5) * dy;
    }

    std::fstream f;
    f.open(_output_directory.path() / "scatterer_distribution.dat", std::ios::out);

    f << "position        population\n";
    for (unsigned i=0; i<pop.size(); ++i){
      f << std::scientific << pos[i] << " " << double(pop[i])/double(_all_scat_list.size()) << "\n";
    }
    f.close();
    
  }

}; // end class monte_carlo

} // end namespace mc

#endif // monte_carlo_h
