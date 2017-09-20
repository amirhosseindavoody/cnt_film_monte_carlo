#include <iostream>
#include <cmath>

#include "ff.h"
#include "particle.h"

namespace mc
{
  // perform free flight within the simulation domain
  void particle::fly(const mc::t_float& dt, const std::pair<mc::arr1d, mc::arr1d>& domain)
  {
    _old_pos = _pos;
    _old_velocity = _velocity;

    _pilot->fly(_pos, _velocity, _eff_mass, dt);
    _pilot->check_boundary(_pos, _velocity, _old_pos, _old_velocity, _eff_mass, dt, domain);
  };

} // namespace mc
