#ifndef utility_h
#define utility_h

#include <cstdlib>
#include <iostream>
#include <list>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <experimental/filesystem>

namespace mc
{

//#################################################################################################
// custom data types
//#################################################################################################
typedef long double t_float; // custom float type for mc class
typedef int t_int;			 // custom integer type for mc class

typedef std::vector<t_float> v1d; // custom 1d float array for mc class
typedef std::vector<v1d> v2d;	 // custom 2d float array for mc class
typedef std::vector<v2d> v3d;	 // custom 3d float array for mc class

typedef std::array<t_float, 3> arr1d; // custom 1d array for expressing the position and momentum in the mc class

//#################################################################################################
// physical and mathematical constants
//#################################################################################################
const t_float eV = 1.6e-19;				  // electron volt in Jouls
const t_float hbar = 6.582e-16 * eV;	  // planck constant in eV.s/rad units
const t_float pi = 3.14159265359;		  // pi number
const t_float elec_mass = 9.10938356e-31; // electron mass
const t_float kB = 1.38064852e-23;		  // Boltzmann constant in Jouls/Kelvin units

//#################################################################################################
// random number routines
//#################################################################################################
void init_random_number_generator(); // initialize the seed of the random number generator

template <typename T>
inline T get_rand_include_zero() // Gets a random number between 0 and 1 including zero.
{
	return static_cast<T>(std::rand()) / static_cast<T>(RAND_MAX);
};

template <typename T>
inline T get_rand_exclude_zero() // Gets a random number between 0 and 1 excluding zero
{
	int r;
	while ((r = rand()) == 0)
	{
	}
	return static_cast<T>(r) / static_cast<T>(RAND_MAX);
};

//#################################################################################################
// helper math functions
//#################################################################################################
template <typename T>
inline T norm(std::array<T, 3> arr) // calculate norm of an array
{
	T result = 0.;
	for (int i = 0; i < arr.size(); ++i)
	{
		result += arr[i] * arr[i];
	}
	return std::sqrt(result);
};

template <typename T>
inline T norm2(std::array<T, 3> arr) // calculate norm^2 of an array
{
	T result = 0.;
	for (int i = 0; i < arr.size(); ++i)
	{
		result += arr[i] * arr[i];
	}
	return result;
};

inline void hist(std::list<mc::t_float> &mlist, const mc::t_int nbin) // UNFINISHED: calculate histogram with nbins
{
	mc::t_float minval = *std::min_element(mlist.begin(), mlist.end());
	mc::t_float maxval = *std::max_element(mlist.begin(), mlist.end());

	std::cout << "minval = " << minval << " , maxval = " << maxval << std::endl;
};

} // end namespace mc

#endif