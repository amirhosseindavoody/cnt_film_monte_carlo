#ifndef forster_particle_h
#define forster_particle_h

#include <iostream>
#include <array>
#include <memory>

#include "particle.h"

namespace mc
{

class forster_particle : public particle
{
private:

protected:

public:
	forster_particle(const mc::arr1d& pos, const mc::t_float& eff_mass,
		const std::shared_ptr<mc::free_flight>& pilot,
		const std::shared_ptr<mc::scatter>& m_scatterer, const mc::t_int& id = -1) // constructor
	{
		set_id(id);
		set_old_pos(pos);
		set_pos(pos);
		set_velocity({0.,0.,0.});
		set_mass(eff_mass);
		set_pilot(pilot);
		set_scatterer(m_scatterer);

		kin_energy();
		_ff_time = scatterer()->ff_time();
	};
	inline void fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // perform free flight within the simulation domain
	{
    pilot()->fly(this, dt);
    pilot()->check_boundary(this, dt, domain);
  };

}; //forster_particle class

} //mc namespace

#endif // forster_particle_h
