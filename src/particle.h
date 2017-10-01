#ifndef particle_h
#define particle_h

#include <iostream>
#include <array>
#include <memory>

#include "scatter.h"
#include "utility.h"

namespace mc
{

class free_flight;
// class scatter;

class particle
{
private:
	mc::t_int _id; // particle id to indentify easily between particles
	mc::arr1d _pos; // position of the particle
	mc::arr1d _velocity; // velocity of the particle
	mc::arr1d _old_pos; // position of the particle in the previous time step, this is used for boundary collision detection
	mc::arr1d _old_velocity; // velocity of the particle in the previous time step, this is used for boundary collision detection
	mc::t_float _eff_mass; // effective mass of the particle

	std::shared_ptr<mc::free_flight> _pilot; // pointer to free_flight object for driving the particle
	std::shared_ptr<mc::scatter> _scatterer; // pointer to scatter object for scattering the particle

protected:
	mc::t_float _ff_time; // free flight time until next scattering event
	mc::t_float _kin_energy; // kinetic energy of the particle

	inline virtual void set_old_pos(const mc::arr1d& value) // set the old_position of the particle
	{
		_old_pos = value;
	}

public:
	inline virtual void set_id(const mc::t_int& id) // set particle id
	{
		_id = id;
	};
	inline virtual const mc::t_int& id() const // get particle id
	{
		return _id;
	};
	inline virtual void set_mass(const mc::t_float value) // set effective mass of the particle
	{
		_eff_mass = value;
	};
	inline virtual const mc::t_float& mass() const // get the effective mass of the particle
	{
		return _eff_mass;
	};
	inline virtual void set_pilot(const std::shared_ptr<mc::free_flight>& pilot) // set the pilot free_flight pointer object
	{
		_pilot = pilot;
	};
	inline virtual const std::shared_ptr<mc::free_flight>& pilot() const // get the pilot free_flight pointer
	{
		return _pilot;
	};
	inline virtual void set_scatterer(const std::shared_ptr<mc::scatter>& scatterer) // set the pointer to the scatterer object
	{
		_scatterer = scatterer;
	};
	inline virtual const std::shared_ptr<mc::scatter>& scatterer() const // return the pointer to the scatterer object
	{
		return _scatterer;
	};
	inline virtual const mc::arr1d& pos() const // get position of the particle
	{
		return _pos;
	};
	inline virtual const mc::t_float& pos(const mc::t_int& i) const // get position of the particle
	{
		return _pos[i];
	};
	inline virtual const mc::arr1d& old_pos() const // get old position of the particle
	{
		return _old_pos;
	};
	inline virtual const mc::t_float& old_pos(const mc::t_int& i) const // get old position of the particle
	{
		return _old_pos[i];
	};
	inline virtual void set_pos(const mc::arr1d& pos) // set position of the particle and set the old position into _old_pos
	{
		_old_pos = _pos;
		_pos = pos;
	};
	inline virtual void set_pos(const mc::t_uint& i, const mc::t_float& value) // set an element of particle position and set the old position into _old_pos
	{
		_old_pos[i] = _pos[i];
		_pos[i] = value;
	};
	inline virtual void rewind_pos() // rewind the current position to the old_pos and do not update the old_pos
	{
		_pos = _old_pos;
	};
	inline virtual const mc::arr1d& velocity() const // get velocity of the particle
	{
		return _velocity;
	};
	inline virtual const mc::t_float& velocity(const mc::t_int& i) const // get velocity of the particle
	{
		return _velocity[i];
	};
	inline virtual const mc::arr1d& old_velocity() const // get the old velocity of the particle
	{
		return _old_velocity;
	};
	inline virtual const mc::t_float& old_velocity(const mc::t_int& i) const // get an element of the old velocity of the particle
	{
		return _old_velocity[i];
	};
	inline virtual void set_velocity(const mc::arr1d& velocity) // set velocity of the particle and put the old velocity into the _old_velocity
	{
		_old_velocity = _velocity;
		_velocity = velocity;
	};
	inline virtual void set_velocity(const mc::t_uint& i, const mc::t_float& value) // set an element of particle velocity and put the old velocity into the _old_velocity
	{
		_old_velocity[i] = _velocity[i];
		_velocity[i] = value;
	};
	inline virtual const mc::t_float& ff_time() const // return the free flight time until the next scattering
	{
		return _ff_time;
	};
	inline virtual mc::t_float& ff_time() // return the free flight time until the next scattering
	{
		return _ff_time;
	};
	inline virtual mc::t_float& get_ff_time() // return the free flight time until the next scattering
	{
		return _ff_time;
	};
	inline virtual void step(mc::t_float dt, const std::pair<mc::arr1d, mc::arr1d>& domain) // step particle state for dt in time
	{
		while(_ff_time <= dt)
		{
			dt -= _ff_time;
			fly(_ff_time, domain);
			_scatterer->update_state(this);
			_ff_time = _scatterer->ff_time();
		}

		fly(dt, domain);
		_ff_time -= dt;
	};
	inline const mc::t_float& kin_energy() // update and return the kinetic energy of the particle
	{
		_kin_energy = mass()*mc::norm2(velocity())/2.;
		return _kin_energy;
	};

	////////////////////////////
	// abstract member functions
	////////////////////////////
	virtual void fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain) = 0; // perform free flight within the simulation domain
	// reinitialize particle properties instead of creating new particles
	virtual void reinitialize(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner, const mc::t_float& beta, const mc::t_float& mass, const std::shared_ptr<mc::free_flight>& pilot, const std::shared_ptr<mc::scatter>& m_scatterer) = 0;

}; //particle class

} //mc namespace

#endif // particle_h
