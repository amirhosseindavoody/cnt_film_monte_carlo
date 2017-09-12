#ifndef scatter_h
#define scatter_h

#include <iostream>
#include <array>

#include "utility.h"

namespace mc
{

class scatter
{
private:
	mc::t_float _max_rate; // maximum scattering rate in the scattering table
	mc::t_float _inverse_max_rate; // inverse of the maximum scattering rate which is the lifetime
	mc::t_int _num_scat_mech; // number of scattering mechanisms
	mc::t_int _num_states; // number of states that the particles can scatter out of.
	mc::t_float _min_energy; // minimum energy for particles considered in the scattering table
	mc::t_float _max_energy; // maximum energy for particles considered in the scattering table
	mc::t_float _delta_energy; // size of energy sections in scattering table
	mc::v2d _scat_table; // scattering table: the first index is the number of states, the second index is the scattering mechanism

	void make_table(); // make the scattering table
	mc::t_int cnvrt_energy_to_state_index(mc::t_float energy); // converts particle energy to appropriate index in the scattering table
public:
	scatter(); // constructor
	inline mc::t_float ff_time() // get random free flight time
	{
		return -_inverse_max_rate*std::log(mc::get_rand_exclude_zero<mc::t_float>());
	};
	inline mc::t_int get_scat_mechanism(mc::t_float energy) // get the index of the scattering mechanism
	{
		mc::t_int state_index = cnvrt_energy_to_state_index(energy);
		mc::t_float random = mc::get_rand_include_zero<mc::t_float>();
		mc::t_int scat_mechanism = 0;
		for (mc::v1d::iterator it=_scat_table[state_index].begin(); it!=_scat_table[state_index].end(); ++it)
		{
			if (*it > random)
			{
				return scat_mechanism;
			}
			scat_mechanism ++;
		}
		return scat_mechanism;
	};
	inline void update_state(const mc::t_int& scat_mechanism, const mc::t_float& energy, mc::arr1d& pos, mc::arr1d& velocity) // update the final state of the particle
	{
		switch(scat_mechanism)
		{
			case 0:
			{
				// std::cout << "0-type scattering" << std::endl;
				break;
			}
			case 1:
			{
				// std::cout << "1-type scattering" << std::endl;
				break;
			}
			case 2:
			{
				// std::cout << "2-type scattering" << std::endl;
				break;
			}
			case 3:
			{
				// std::cout << "3-type scattering" << std::endl;
				break;
			}
			case 4:
			{
				// std::cout << "4-type scattering" << std::endl;
				break;
			}
			case 5:
			{
				// std::cout << "5-type scattering" << std::endl;
				break;
			}
			case 6:
			{
				// std::cout << "6-type scattering" << std::endl;
				break;
			}
			case 7:
			{
				// std::cout << "7-type scattering" << std::endl;
				break;
			}
			case 8:
			{
				// std::cout << "8-type scattering" << std::endl;
				break;
			}
			case 9:
			{
				// std::cout << "9-type scattering" << std::endl;
				break;
			}
		}

		// get random energy with correct distribution
		const mc::t_float max_energy = 10.*0.025*mc::eV;
		mc::t_float new_energy = energy;
		if (new_energy > max_energy) // renormalize energy if the electron has too much energy
		{
			new_energy = -(3./2.*0.025*mc::eV)*std::log(mc::get_rand_include_zero<mc::t_float>());
		}
		mc::t_float velocity_magnitude = std::sqrt(new_energy*2./mc::elec_mass);


		// get uniformly distribution direction
		mc::t_float theta = std::acos(1.-2.*mc::get_rand_include_zero<mc::t_float>());
		mc::t_float phi = 2.*mc::pi*mc::get_rand_include_zero<mc::t_float>();
		velocity = {velocity_magnitude*std::sin(theta)*std::cos(phi), velocity_magnitude*std::sin(theta)*std::sin(phi), velocity_magnitude*std::cos(theta)};

		// pos = {0., 0., 0.};
	};

}; // end class scatter

} // end namespace mc

#endif // scatter_h
