#include <iostream>
#include <array>
#include <armadillo>

#include "../helper/utility.h"
#include "discrete_forster_ff.h"

namespace mc
{

void discrete_forster_free_flight::check_boundary(discrete_forster_free_flight::t_particle* p, const double &dt, const std::pair<arma::vec, arma::vec>& domain) // check for collision to boundaries
{
	if (arma::any(p->pos()<domain.first)) {
		p->rewind_pos();
		return;
	}

	if (arma::any(p->pos()>domain.second)) {
		p->rewind_pos();
		return;
	}
};

} // end namespace mc
