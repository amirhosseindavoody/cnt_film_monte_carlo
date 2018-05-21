#ifndef scatterer_h
#define scatterer_h

#include <iostream>
#include <array>
#include <fstream>
#include <experimental/filesystem>
#include <regex>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <tuple>
#include <armadillo>

#include "../helper/utility.h"

#include "./discrete_forster_particle.h"
#include "./scattering_struct.h"

namespace mc
{

class scatterer
{
  struct neighbor {
    neighbor(std::shared_ptr<scatterer> s, double d,
             double rate)
        : s_ptr(s), distance(d), rate(d){};
    std::shared_ptr<scatterer> s_ptr;
    double distance;
    double rate;
  };

private:
	double _max_rate; // maximum scattering rate in the scattering table
	double _inverse_max_rate; // inverse of the maximum scattering rate which is the lifetime
	arma::vec _pos; // position of the scatterer
  arma::vec _orientation; // orientation of the scattering object site
  std::list<scatterer::neighbor> _neighbors;

public:
  std::vector< std::vector<scatterer*>* > close_scats; // index of neighboring grid cells
  const scattering_struct* scat_tab;

	// default constructor
  scatterer(){};

  // set position of the scatterer
  void set_pos(const arma::vec& position) {
    _pos = position;
  };

  // set a component of the scatterer position
  void set_pos(const unsigned& i, const double& value) {
    _pos(i) = value;
  };

  // get position of the scatterer
  const arma::vec& pos() const {
    return _pos;
  };

  // get i'th component of position of the scatterer
  const double& pos(const unsigned& i) const { return _pos(i); };

  // set the orientation of the scatterer object
  void set_orientation(const arma::vec& m_orientation) {
    _orientation = m_orientation;
  };

  // set the i'th element of the orientation of the scatterer object
  void set_orientation(const mc::t_uint& i, const double& value) {
    _orientation(i) = value;
  };

  // get the orientation of the scatterer object
  const arma::vec& orientation() const { return _orientation; };

  // get the i'th element of the orientation of the scatterer object
  const double& orientation(const unsigned& i) const {
    return _orientation(i);
  };

  // get random free flight time
  double ff_time() const {
    int r;
    while ((r = rand()) == 0) {
    }

    return -_inverse_max_rate * std::log(double(r) / double(RAND_MAX));
  };

  // update the final state of the particle
  void update_state(discrete_forster_particle* p,
                    const double& max_hopping_radius);

  // add a scattering object and its distance to the neighbors list
  void add_neighbor(const std::shared_ptr<scatterer>& neighbor_ptr,
                    const double& distance, const double& rate) {
    _neighbors.emplace_back(neighbor_ptr, distance, rate);
  };

  // get the number of neighbors
  unsigned number_of_neighbors() { return _neighbors.size(); };

  // sort neighbors in the ascending order of distances
  void sort_neighbors() {
    _neighbors.sort([](const auto& n1, const auto& n2) {
      return n1.distance < n2.distance;
    });
  };

  // make a cumulative scattering rate table
  void make_cumulative_scat_rate() {
		_max_rate = 0;
		for (auto& neighbor: _neighbors)
		{
      _max_rate += neighbor.rate;
      neighbor.rate = _max_rate;
		}
		_inverse_max_rate = 1/_max_rate;
	};

  // print neighbor distances
  void print_neighbor_distances() {
		std::cout << "...neighbor distances: ";
		for (const auto& neighbor : _neighbors)
		{
      std::cout << neighbor.distance << " , ";
		}
		std::cout << std::endl;
	};

  // print neighbor rates
  void print_neighbor_rates() {
		std::cout << "...neighbor rates: ";
		for (const auto& neighbor : _neighbors)
		{
      std::cout << neighbor.rate << " , ";
		}
		std::cout << std::endl;
	};

  // return closest neighbor distance
  double closest_neighbor() {
    auto cmp_neighbor = [](const scatterer::neighbor& n1,
                           const scatterer::neighbor& n2) {
      return n1.distance < n2.distance;
    };
    return std::min_element(_neighbors.begin(),_neighbors.end(),cmp_neighbor)->distance;
  };

  // return farthest neighbor distance
  double farthest_neighbor() {
    auto cmp_neighbor = [](const scatterer::neighbor& n1,
                           const scatterer::neighbor& n2) {
      return n1.distance < n2.distance;
    };
    return std::max_element(_neighbors.begin(),_neighbors.end(),cmp_neighbor)->distance;
  };

  // return the maximum scattering rate
  const double& max_rate() const { return _max_rate; };

 private:
  // find neighbors of the current scatterer and their scattering rates
  std::vector < std::pair<double, scatterer*>> find_neighbors(const double& max_hopping_radius);

};  // end class scatterer

} // end namespace mc

#endif  // scatterer_h