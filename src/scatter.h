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
	mc::t_float ff_time(); // get random free flight time
	mc::t_int get_scat_mechanism(mc::t_float energy); // get the index of the scattering mechanism
	void update_state(const mc::t_int& scat_mechanism, const mc::t_float& energy, mc::arr1d& pos, mc::arr1d& velocity); // update the final state of the particle

}; // end class scatter

} // end namespace mc

#endif // scatter_h