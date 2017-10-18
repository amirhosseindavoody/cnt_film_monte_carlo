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
public:
  typedef mc::discrete_forster_particle t_particle; // particle type

private:
	mc::t_float _max_rate; // maximum scattering rate in the scattering table
	mc::t_float _inverse_max_rate; // inverse of the maximum scattering rate which is the lifetime
	mc::arr1d _pos; // position of the scatterer
	std::list<std::tuple<std::shared_ptr<mc::discrete_forster_scatter>, mc::t_float, mc::t_float>> _neighbors; // list of pairs containing distance and shared pointer to the neighboring scatterers

	// compare components of neighbors list by measure of their distance
	static bool _compare_neighbors(const std::tuple<std::shared_ptr<mc::discrete_forster_scatter>, mc::t_float, mc::t_float>& lhs, const std::tuple<std::shared_ptr<mc::discrete_forster_scatter>, mc::t_float, mc::t_float>& rhs)
	{
		return std::get<1>(lhs) < std::get<1>(rhs);
	};

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
		_neighbors.push_back(std::make_tuple(neighbor_ptr, distance, 1.e12*std::pow(1.e-9,6)/std::pow(distance,6)));
	};
  // get the number of neighbors
  mc::t_uint number_of_neighbors()
	{
		return _neighbors.size();
	};
  // sort neighbors in the ascending order
  void sort_neighbors()
	{
		_neighbors.sort(_compare_neighbors);
	};
  // make a cumulative scattering rate table
  void make_cumulative_scat_rate()
	{
		_max_rate = 0;
		for (auto& neighbor: _neighbors)
		{
			_max_rate += std::get<2>(neighbor);
			std::get<2>(neighbor) = _max_rate;

			// std::cout << "distance = " << std::get<1>(neighbor) << "....rate = " << std::get<2>(neighbor) << std::endl;
			// std::cin.ignore();
		}
		_inverse_max_rate = 1/_max_rate;

		// std::cout << "number of neighbors = " << _neighbors.size() << std::endl;
		// std::cout << "max scat rate: " << _max_rate << std::endl;
		// std::cin.ignore();
	};
  // print neighbor distances
  void print_neighbor_distances()
	{
		std::cout << "...neighbor distances: ";
		for (const auto& neighbor : _neighbors)
		{
			std::cout << std::get<2>(neighbor) << " , ";
		}
		std::cout << std::endl;
	};

}; // end class discrete_forster_scatter


} // end namespace mc

#endif // discrete_forster_scatter_h
