#ifndef discrete_forster_particle_h
#define discrete_forster_particle_h

#include <iostream>
#include <array>
#include <memory>
#include <armadillo>

#include "discrete_forster_ff.h"
#include "scatterer.h"

namespace mc
{

class discrete_forster_particle
{
public:
	typedef mc::discrete_forster_particle t_particle; // particle type
	typedef mc::discrete_forster_free_flight t_ff; // free_flight type
	// typedef mc::discrete_forster_region t_region; // region type

private:
	arma::vec _pos; // position of the particle
	arma::vec _old_pos; // position of the particle in the previous time step, this is used for boundary collision detection

	std::shared_ptr<t_ff> _pilot; // pointer to free_flight object for driving the particle
	std::shared_ptr<scatterer> _scat_ptr; // pointer to scatter object for scattering the particle

	double _ff_time; // free flight time until next scattering event

public:
	//constructor
  discrete_forster_particle(
      const arma::vec& pos, const std::shared_ptr<t_ff>& pilot,
      const std::shared_ptr<scatterer>& m_scatterer);  // constructor

  // reinitialize particle properties instead of creating new particles
  void reinitialize(const arma::vec& lower_corner, const arma::vec& upper_corner,
                    const double& beta, const double& mass,
                    const std::shared_ptr<t_ff>& pilot,
                    const std::shared_ptr<scatterer>& m_scatterer);

  // perform free flight within the simulation domain
  void fly(const double& dt,
           const std::pair<arma::vec, arma::vec>&
               domain);  // perform free flight within the simulation domain

  // set the pilot free_flight pointer object
  void set_pilot(const std::shared_ptr<t_ff>& pilot) {
    _pilot = pilot;
  };

  // get the pilot free_flight pointer
  const std::shared_ptr<t_ff>& pilot() const {
    return _pilot;
  };

  // set the pointer to the scatterer object
  void set_scatterer(const std::shared_ptr<scatterer>& s) {
    _scat_ptr = s;
	};

	// return the pointer to the scatterer object
	const std::shared_ptr<scatterer>& scat_ptr() const {
		return _scat_ptr;
	};

	// get position of the particle
	const arma::vec& pos() const  {
		return _pos;
	};

	// get position of the particle
	const double& pos(const double& i) const  {
		return _pos(i);
	};

	// get old position of the particle
	const arma::vec& old_pos() const  {
		return _old_pos;
	};

	// get old position of the particle
	const double& old_pos(const int& i) const  {
		return _old_pos(i);
	};

	// set position of the particle and set the old position into _old_pos
	void set_pos(const arma::vec& pos) {
		_old_pos = _pos;
		_pos = pos;
	};

	// set an element of particle position and set the old position into _old_pos
	void set_pos(const int& i, const double& value)  {
		_old_pos(i) = _pos(i);
		_pos(i) = value;
	};
	
	// rewind the current position to the old_pos and do not update the old_pos
	void rewind_pos() {
		_pos = _old_pos;
	};

	// set the old position of the particle.
	void set_old_pos(const arma::vec& old_pos) {
		_old_pos = old_pos;
	};

	// return the free flight time until the next scattering
	const double& ff_time() const  {
		return _ff_time;
	};

	// return the free flight time until the next scattering
	void set_ff_time(const double& value)  {
		_ff_time = value;
	};

	// update the _ff_time by calling the underlying scatterer
	void get_ff_time();

	// step particle state for dt in time
	void step(double dt, const std::pair<arma::vec, arma::vec>& domain);

}; //discrete_forster_particle class

} //mc namespace

#endif // discrete_forster_particle_h
