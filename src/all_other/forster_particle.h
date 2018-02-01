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
	// specialized constructor
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

	// a common form constructor
	forster_particle(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner, const mc::t_float& beta, const mc::t_float& mass,
		const std::shared_ptr<mc::free_flight>& pilot,
		const std::shared_ptr<mc::scatter>& m_scatterer) // constructor
	{
		// set_id(id);

		for (unsigned i=0; i<lower_corner.size(); ++i)
		{
			set_pos(i, lower_corner[i]+(upper_corner[i]-lower_corner[i])*mc::get_rand_include_zero<mc::t_float>());
		}
		// set_velocity(mc::rand_velocity(beta, mc::elec_mass));
		set_mass(mass);
		set_pilot(pilot);
		set_scatterer(m_scatterer);

		// kin_energy();
		_ff_time = scatterer()->ff_time();
	};

	// reinitialize particle properties instead of creating new particles
	void reinitialize(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner, const mc::t_float& beta, const mc::t_float& mass, const std::shared_ptr<mc::free_flight>& pilot, const std::shared_ptr<mc::scatter>& m_scatterer)
	{
		for (unsigned i=0; i<lower_corner.size(); ++i)
		{
			set_pos(i, lower_corner[i]+(upper_corner[i]-lower_corner[i])*mc::get_rand_include_zero<mc::t_float>());
		}
		_ff_time = scatterer()->ff_time();
	};

	// perform free flight within the simulation domain
	inline void fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // perform free flight within the simulation domain
	{
    // pilot()->fly(this, dt);
    pilot()->check_boundary(this, dt, domain);
  };

}; //forster_particle class

} //mc namespace

#endif // forster_particle_h
