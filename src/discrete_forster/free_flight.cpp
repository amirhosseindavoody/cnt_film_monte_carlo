#include <iostream>
#include <array>
#include <armadillo>

#include "../helper/utility.h"
#include "free_flight.h"

namespace mc
{

// check for collision to boundaries
void free_flight::check_boundary(particle* p, const double& dt, const std::pair<arma::vec, arma::vec>& domain) const {
  if (arma::any(p->pos() < domain.first)) {
    p->rewind_pos();
    return;
  }

  if (arma::any(p->pos() > domain.second)) {
    p->rewind_pos();
    return;
  }
};

} // end namespace mc
