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
	// a specialized constructor
	gas_particle(const mc::arr1d& pos, const mc::arr1d& velocity, const mc::t_float& eff_mass,
		const std::shared_ptr<mc::free_flight>& pilot,
		const std::shared_ptr<mc::scatter>& m_scatterer, const mc::t_int& id = -1) // constructor
	{
		set_id(id);
		set_pos(pos);
		set_velocity(velocity);
		set_mass(eff_mass);
		set_pilot(pilot);
		set_scatterer(m_scatterer);

		kin_energy();
		_ff_time = scatterer()->ff_time();
	};

	// a common form constructor
	gas_particle(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner, const mc::t_float& beta, const mc::t_float& mass,
		const std::shared_ptr<mc::free_flight>& pilot,
		const std::shared_ptr<mc::scatter>& m_scatterer) // constructor
	{
		// set_id(id);

		for (unsigned i=0; i<lower_corner.size(); ++i)
		{
			set_pos(i, lower_corner[i]+(upper_corner[i]-lower_corner[i])*double(std::rand())/double(RAND_MAX));
		}
		set_velocity(mc::rand_velocity(beta, mc::elec_mass));
		set_mass(mass);
		set_pilot(pilot);
		set_scatterer(m_scatterer);

		kin_energy();
		_ff_time = scatterer()->ff_time();
	};

	// reinitialize particle properties instead of creating new particles
	void reinitialize(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner, const mc::t_float& beta, const mc::t_float& mass, const std::shared_ptr<mc::free_flight>& pilot, const std::shared_ptr<mc::scatter>& m_scatterer)
	{
		for (unsigned i=0; i<lower_corner.size(); ++i)
		{
			set_pos(i, lower_corner[i]+(upper_corner[i]-lower_corner[i])*double(std::rand())/double(RAND_MAX));
		}
		set_velocity(mc::rand_velocity(beta, mc::elec_mass));
		kin_energy();
		_ff_time = scatterer()->ff_time();
	};

	// perform free flight within the simulation domain
	inline void fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // perform free flight within the simulation domain
	{
    pilot()->fly(this, dt);
    pilot()->check_boundary(this, dt, domain);
  };

}; //gas_particle class

} //mc namespace

#endif // gas_particle_h
