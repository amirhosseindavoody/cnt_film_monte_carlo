#ifndef gas_particle_h
#define gas_particle_h

#include <iostream>
#include <array>
#include <memory>

#include "particle.h"

namespace mc
{

class gas_particle : public particle
{
private:

protected:

public:
	gas_particle(const mc::arr1d& pos, const mc::arr1d& velocity, const mc::t_float& eff_mass,
		const std::shared_ptr<mc::free_flight>& pilot,
		const std::shared_ptr<mc::scatter>& scatterer, const mc::t_int& id = -1) // constructor
	{
		set_id(id);
		set_pos(pos);
		set_velocity(velocity);
		set_mass(eff_mass);
		set_pilot(pilot);
		set_scatterer(scatterer);

		kin_energy();
		update_ff_time();
	};
	inline void fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // perform free flight within the simulation domain
	{
    pilot()->fly(this, dt);
    pilot()->check_boundary(this, dt, domain);
  };

}; //gas_particle class

} //mc namespace

#endif // gas_particle_h
