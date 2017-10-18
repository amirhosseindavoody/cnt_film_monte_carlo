#ifndef discrete_forster_particle_h
#define discrete_forster_particle_h

#include <iostream>
#include <array>
#include <memory>

// #include "all_particles.h"
#include "discrete_forster_scatter.h"
#include "discrete_forster_ff.h"

namespace mc
{

class discrete_forster_particle
{
public:
	typedef mc::discrete_forster_particle t_particle; // particle type
	typedef mc::discrete_forster_free_flight t_ff; // free_flight type
	// typedef mc::discrete_forster_region t_region; // region type
	typedef mc::discrete_forster_scatter t_scatter; // scatter

private:
	mc::arr1d _pos; // position of the particle
	mc::arr1d _old_pos; // position of the particle in the previous time step, this is used for boundary collision detection

	std::shared_ptr<t_ff> _pilot; // pointer to free_flight object for driving the particle
	std::shared_ptr<t_scatter> _scatterer; // pointer to scatter object for scattering the particle

	mc::t_float _ff_time; // free flight time until next scattering event

public:
	//constructor
	discrete_forster_particle(const mc::arr1d& pos, const std::shared_ptr<t_ff>& pilot, const std::shared_ptr<t_scatter>& m_scatterer); // constructor

	// reinitialize particle properties instead of creating new particles
	void reinitialize(const mc::arr1d& lower_corner, const mc::arr1d& upper_corner, const mc::t_float& beta, const mc::t_float& mass, const std::shared_ptr<t_ff>& pilot, const std::shared_ptr<t_scatter>& m_scatterer);

	// perform free flight within the simulation domain
	void fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain); // perform free flight within the simulation domain

	void set_pilot(const std::shared_ptr<t_ff>& pilot) // set the pilot free_flight pointer object
	{
		_pilot = pilot;
	};
	const std::shared_ptr<t_ff>& pilot() const // get the pilot free_flight pointer
	{
		return _pilot;
	};
	void set_scatterer(const std::shared_ptr<t_scatter>& scatterer) // set the pointer to the scatterer object
	{
		_scatterer = scatterer;
	};
	const std::shared_ptr<t_scatter>& scatterer() const // return the pointer to the scatterer object
	{
		return _scatterer;
	};
	const mc::arr1d& pos() const // get position of the particle
	{
		return _pos;
	};
	const mc::t_float& pos(const mc::t_int& i) const // get position of the particle
	{
		return _pos[i];
	};
	const mc::arr1d& old_pos() const // get old position of the particle
	{
		return _old_pos;
	};
	const mc::t_float& old_pos(const mc::t_int& i) const // get old position of the particle
	{
		return _old_pos[i];
	};
	void set_pos(const mc::arr1d& pos) // set position of the particle and set the old position into _old_pos
	{
		_old_pos = _pos;
		_pos = pos;
	};
	void set_pos(const mc::t_uint& i, const mc::t_float& value) // set an element of particle position and set the old position into _old_pos
	{
		_old_pos[i] = _pos[i];
		_pos[i] = value;
	};
	void rewind_pos() // rewind the current position to the old_pos and do not update the old_pos
	{
		_pos = _old_pos;
	};
	// set the old position of the particle.
	void set_old_pos(const mc::arr1d& old_pos)
	{
		_old_pos = old_pos;
	};
	const mc::t_float& ff_time() const // return the free flight time until the next scattering
	{
		return _ff_time;
	};
	void set_ff_time(const mc::t_float& value) // return the free flight time until the next scattering
	{
		_ff_time = value;
	};
	// update the _ff_time by calling the underlying scatterer
	void get_ff_time();
	// step particle state for dt in time
	void step(mc::t_float dt, const std::pair<mc::arr1d, mc::arr1d>& domain);


}; //discrete_forster_particle class

} //mc namespace

#endif // discrete_forster_particle_h
