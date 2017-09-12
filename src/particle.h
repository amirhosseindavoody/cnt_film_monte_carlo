#ifndef particle_h
#define particle_h

#include <iostream>
#include <array>
#include <memory>

#include "ff.h"
#include "scatter.h"
#include "utility.h"

namespace mc
{

class particle
{
private:
	mc::t_int _id; // particle id to indentify easily between particles
	mc::arr1d _pos; // position of the particle
	mc::arr1d _velocity; // velocity of the particle
	mc::arr1d _old_pos; // position of the particle in the previous time step, this is used for boundary collision detection
	mc::arr1d _old_velocity; // velocity of the particle in the previous time step, this is used for boundary collision detection
	mc::t_float _eff_mass; // effective mass of the particle
	mc::t_float _kin_energy; // kinetic energy of the particle
	mc::t_float _ff_time; // free flight time until next scattering event

	std::shared_ptr<mc::free_flight> _pilot; // pointer to free_flight object for driving the particle
	std::shared_ptr<mc::scatter> _scatterer; // pointer to scatter object for scattering the particle

public:
	inline particle(const mc::arr1d& pos={0.,0.,0.}, const mc::arr1d& velocity={0.,0.,0.}, const mc::t_float& eff_mass = mc::elec_mass,
		const std::shared_ptr<mc::free_flight>& pilot = std::make_shared<mc::free_flight>(),
		const std::shared_ptr<mc::scatter>& scatterer = std::make_shared<mc::scatter>(), const mc::t_int& id = -1) // constructor
	{
		_id = id;
		_pos = pos;
		_velocity = velocity;
		_eff_mass = eff_mass;
		kin_energy();
		_pilot = pilot;
		_scatterer = scatterer;
		update_ff_time();
	};
	inline void fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // perform free flight within the simulation domain
	{
		_old_pos = _pos;
		_old_velocity = _velocity;

		_pilot->fly(_pos, _velocity, _eff_mass, dt);
		_pilot->check_boundary(_pos, _velocity, _old_pos, _old_velocity, _eff_mass, dt, domain);
	};
	inline const mc::arr1d& pos() const // get position of the particle
	{
		return _pos;
	};
	inline const mc::t_float& pos(const mc::t_int& i) const // get position of the particle
	{
		return _pos[i];
	};
	inline void set_pos(const mc::arr1d& pos) // set position of the particle
	{
		_pos = pos;
	};
	inline const mc::arr1d& velocity() const // get velocity of the particle
	{
		return _velocity;
	};
	inline const mc::t_float& velocity(const mc::t_int& i) const // get velocity of the particle
	{
		return _velocity[i];
	};
	inline void set_velocity(const mc::arr1d& velocity) // set velocity of the particle
	{
		_velocity = velocity;
	};
	inline const mc::t_float& kin_energy() // update and return the kinetic energy of the particle
	{
		_kin_energy = _eff_mass*mc::norm2(_velocity)/2.;
		return _kin_energy;
	};
	inline const mc::t_float& update_ff_time() // update and return the free flight time until the next scattering
	{
		_ff_time = _scatterer->ff_time();
		return _ff_time;
	};
	inline mc::t_float& get_ff_time() // return the free flight time until the next scattering
	{
		return _ff_time;
	};
	inline void scatter() // scatter the particle to a new state using the scatterer object
	{
		mc::t_int scat_mechanism = _scatterer->get_scat_mechanism(kin_energy());
		_scatterer->update_state(scat_mechanism, kin_energy(), _pos, _velocity);
	};
	inline friend std::ostream& operator<< (std::ostream& stream, const particle& _particle) // print the state of the particle
	{
		for (int i=0; i<_particle._pos.size(); ++i)
		{
			stream << _particle._pos[i] << " ";
		}
		stream << ", ";
		for (int i=0; i<_particle._velocity.size(); ++i)
		{
			stream << _particle._velocity[i] << " ";
		}
	};
	inline const mc::t_int& id() const // get a constant reference to the particle id;
	{
		return _id;
	};
	inline void step(mc::t_float dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // step particle state for dt in time
	{
		while(_ff_time <= dt)
		{
			dt -= _ff_time;
			fly(_ff_time, domain);
			scatter();
			update_ff_time();
		}

		fly(dt, domain);
		_ff_time -= dt;
	};

}; //particle class

} //mc namespace

#endif // particle_h
