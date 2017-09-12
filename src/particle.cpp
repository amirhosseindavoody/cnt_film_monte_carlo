#include <iostream>
#include <cmath>

#include "particle.h"

namespace mc
{

// // constructor
// particle::particle(const mc::arr1d pos, const mc::arr1d velocity, const mc::t_float eff_mass, const std::shared_ptr<mc::free_flight> pilot, const std::shared_ptr<mc::scatter> scatterer, const mc::t_int id)
// {
// 	_id = id;
// 	_pos = pos;
// 	_velocity = velocity;
// 	_eff_mass = eff_mass;
// 	kin_energy();
// 	_pilot = pilot;
// 	_scatterer = scatterer;
// 	update_ff_time();
// };

// // perform free flight
// void particle::fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain)
// {
// 	_old_pos = _pos;
// 	_old_velocity = _velocity;
//
// 	_pilot->fly(_pos, _velocity, _eff_mass, dt);
// 	_pilot->check_boundary(_pos, _velocity, _old_pos, _old_velocity, _eff_mass, dt, domain);
// 	// std::cin.ignore();
// };

// // print the state of the particle
// std::ostream& operator<< (std::ostream& stream, const particle& _particle)
// {
// 	for (int i=0; i<_particle._pos.size(); ++i)
// 	{
// 		stream << _particle._pos[i] << " ";
// 	}
// 	stream << ", ";
// 	for (int i=0; i<_particle._velocity.size(); ++i)
// 	{
// 		stream << _particle._velocity[i] << " ";
// 	}
// };

// // update and return the kinetic energy of the particle
// const mc::t_float& particle::kin_energy()
// {
// 	_kin_energy = _eff_mass*mc::norm2(_velocity)/2.;
// 	return _kin_energy;
// };

// // update and return the free flight time until the next scattering
// const mc::t_float& particle::update_ff_time()
// {
// 	_ff_time = _scatterer->ff_time();
// 	return _ff_time;
// };

// // return the free flight time until the next scattering
// mc::t_float& particle::get_ff_time()
// {
// 	return _ff_time;
// };

// // scatter the particle to a new state using the scatterer object
// void particle::scatter()
// {
// 	mc::t_int scat_mechanism = _scatterer->get_scat_mechanism(kin_energy());
// 	_scatterer->update_state(scat_mechanism, kin_energy(), _pos, _velocity);
// };

// // get a constant reference to the particle id;
// const mc::t_int& particle::id() const
// {
// 	return _id;
// };

// // step particle state for dt in time
// void particle::step(mc::t_float dt, const std::pair<mc::arr1d, mc::arr1d>& domain)
// {
// 	while(_ff_time <= dt)
// 	{
// 		dt -= _ff_time;
// 		fly(_ff_time, domain);
// 		scatter();
// 		update_ff_time();
// 	}
//
// 	fly(dt, domain);
// 	_ff_time -= dt;
// };

} // namespace mc
