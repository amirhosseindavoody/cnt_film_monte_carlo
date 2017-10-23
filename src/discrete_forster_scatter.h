#ifndef discrete_forster_scatter_h
#define discrete_forster_scatter_h

#include <iostream>
#include <array>
#include <fstream>
#include <experimental/filesystem>
#include <regex>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <tuple>

#include "utility.h"
#include "discrete_forster_particle.h"

namespace mc
{

class discrete_forster_scatter
{
  struct neighbor
  {
    neighbor(std::shared_ptr<mc::discrete_forster_scatter> m_scatterer, mc::t_float m_distance, mc::t_float m_rate): scatterer(m_scatterer), distance(m_distance), rate(m_rate)
    {};
    std::shared_ptr<mc::discrete_forster_scatter> scatterer;
    mc::t_float distance;
    mc::t_float rate;
    // bool operator <(const neighbor& other)
    // {
    //   return distance < other.distance;
    // };
  };

public:
  typedef mc::discrete_forster_particle t_particle; // particle type

private:
	mc::t_float _max_rate; // maximum scattering rate in the scattering table
	mc::t_float _inverse_max_rate; // inverse of the maximum scattering rate which is the lifetime
	mc::arr1d _pos; // position of the scatterer
  mc::arr1d _orientation; // orientation of the scattering object site
  std::list<mc::discrete_forster_scatter::neighbor> _neighbors;

public:

	// constructor
	discrete_forster_scatter()
  {};
  // set position of the scatterer
  void set_pos(const mc::arr1d& position)
	{
		_pos = position;
	};
  // set a component of the scatterer position
  void set_pos(const mc::t_uint& i, const mc::t_float& value)
	{
		_pos[i] = value;
	}
  // get position of the scatterer
  const mc::arr1d& pos() const
	{
		return _pos;
	};
  // get i'th component of position of the scatterer
  const mc::t_float& pos(const mc::t_uint& i) const
	{
		return _pos[i];
	};
  // set the orientation of the scatterer object
  void set_orientation(const mc::arr1d& m_orientation)
  {
    _orientation = m_orientation;
  };
  // set the i'th element of the orientation of the scatterer object
  void set_orientation(const mc::t_uint& i, const mc::t_float& value)
  {
    _orientation[i] = value;
  };
  // get the orientation of the scatterer object
  const mc::arr1d& orientation() const
  {
    return _orientation;
  };
  // get the i'th element of the orientation of the scatterer object
  const mc::t_float& orientation(const mc::t_uint& i) const
  {
    return _orientation[i];
  };
  // get random free flight time
  mc::t_float ff_time() const
	{
		return -_inverse_max_rate*std::log(mc::get_rand_exclude_zero<mc::t_float>());
	};
	// update the final state of the particle
	void update_state(t_particle* p);
	// add a scattering object and its distance to the neighbors list
	void add_neighbor(const std::shared_ptr<mc::discrete_forster_scatter>& neighbor_ptr, const mc::t_float& distance)
	{
    mc::t_float rate = 1.e12*std::pow(1.e-9,6)/std::pow(distance,6);
    _neighbors.emplace_back(neighbor_ptr, distance, rate);
	};
  // get the number of neighbors
  mc::t_uint number_of_neighbors()
	{
		return _neighbors.size();
	};
  // sort neighbors in the ascending order
  void sort_neighbors()
	{
    auto cmp_neighbor = [](const mc::discrete_forster_scatter::neighbor& n1, const mc::discrete_forster_scatter::neighbor& n2)
    {
      return n1.distance < n2.distance;
    };
    _neighbors.sort(cmp_neighbor);
	};
  // make a cumulative scattering rate table
  void make_cumulative_scat_rate()
	{
		_max_rate = 0;
		for (auto& neighbor: _neighbors)
		{
      _max_rate += neighbor.rate;
      neighbor.rate = _max_rate;
		}
		_inverse_max_rate = 1/_max_rate;
	};
  // print neighbor distances
  void print_neighbor_distances()
	{
		std::cout << "...neighbor distances: ";
		for (const auto& neighbor : _neighbors)
		{
      std::cout << neighbor.distance << " , ";
		}
		std::cout << std::endl;
	};
  // print neighbor rates
  void print_neighbor_rates()
	{
		std::cout << "...neighbor rates: ";
		for (const auto& neighbor : _neighbors)
		{
      std::cout << neighbor.rate << " , ";
		}
		std::cout << std::endl;
	};
  // return closest neighbor distance
  mc::t_float closest_neighbor()
  {
    auto cmp_neighbor = [](const mc::discrete_forster_scatter::neighbor& n1, const mc::discrete_forster_scatter::neighbor& n2)
    {
      return n1.distance < n2.distance;
    };
    return std::min_element(_neighbors.begin(),_neighbors.end(),cmp_neighbor)->distance;
  };
  // return farthest neighbor distance
  mc::t_float farthest_neighbor()
  {
    auto cmp_neighbor = [](const mc::discrete_forster_scatter::neighbor& n1, const mc::discrete_forster_scatter::neighbor& n2)
    {
      return n1.distance < n2.distance;
    };
    return std::max_element(_neighbors.begin(),_neighbors.end(),cmp_neighbor)->distance;
  };
  // return the maximum scattering rate
  const mc::t_float& max_rate() const
  {
    return _max_rate;
  };

}; // end class discrete_forster_scatter


} // end namespace mc

#endif // discrete_forster_scatter_h
