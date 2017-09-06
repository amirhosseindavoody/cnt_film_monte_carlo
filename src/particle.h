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
	particle(const mc::arr1d _pos={0.,0.,0.}, const mc::arr1d _velocity={0.,0.,0.}, const mc::t_float _eff_mass = mc::elec_mass, const std::shared_ptr<mc::free_flight> _pilot = std::make_shared<mc::free_flight>(), const std::shared_ptr<mc::scatter> _scatterer = std::make_shared<mc::scatter>(), const mc::t_int id = -1); // constructor
	void fly(const mc::t_float& dt, const mc::arr1d& volume); // perform free flight
	const mc::arr1d& pos(); // get position of the particle
	const mc::t_float& pos(mc::t_int i); // get position of the particle
	const mc::arr1d& velocity(); // get velocity of the particle
	const mc::t_float& kin_energy(); // update and return the kinetic energy of the particle
	const mc::t_float& update_ff_time(); // update and return the free flight time until the next scattering
	mc::t_float& get_ff_time(); // return the free flight time until the next scattering
	void scatter(); // scatter the particle to a new state using the scatterer object
	friend std::ostream& operator<< (std::ostream& stream, const particle& _particle); // print the state of the particle
	const mc::t_int& id() const; // get a constant reference to the particle id;

}; //particle class

} //mc namespace

#endif // particle_h
