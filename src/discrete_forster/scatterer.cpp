#include <iostream>
#include <limits>

#include "discrete_forster_particle.h"
#include "scatterer.h"

namespace mc
{
  // update the final state of the particle
  void scatterer::update_state(discrete_forster_particle* p, const double& max_hopping_radius) {

    using namespace std;
    vector<pair<double, scatterer*>> nlist = find_neighbors(max_hopping_radius);

    if (!nlist.empty()){
      double dice = nlist.back().first * double(std::rand()) / double(RAND_MAX);
      unsigned left = -1;
      unsigned right = nlist.size()-1;

      while(left+1<right){
        unsigned mid = (left+right)/2;
        if (nlist[mid].first<=dice){
          left=mid;
        } else {
          right=mid;
        }
      }

      p->set_pos(nlist[right].second->pos());
      p->set_scatterer(nlist[right].second);
    }
	};

  // find neighbors of the current scatterer and their scattering rates
  std::vector<std::pair<double, scatterer*>>
  scatterer::find_neighbors(const double& max_hopping_radius){

    std::vector<std::pair<double, scatterer*>> neighbors_list;

    auto check_neighbor = [&max_hopping_radius, this, &neighbors_list]
                          (scatterer& s1, scatterer& s2) {
      double cosTheta, theta, y1, y2, sin2Theta, axis_shift_1, axis_shift_2,
          z_shift;

      arma::vec dR = s1.pos() - s2.pos();
      double distance = arma::norm(dR);

      if ((distance < max_hopping_radius) && (distance > 0.4e-9)) {
        arma::vec a1 = s1.orientation();
        arma::vec a2 = s2.orientation();
        cosTheta = arma::dot(a1, a2);
        // check if parallel case has happend
        if (cosTheta == 1) {
          axis_shift_1 = 0;
          axis_shift_2 = arma::dot(dR, a1);
          theta = 0;
          z_shift = arma::norm(dR - arma::dot(dR, a1) * a1);
        } else {
          theta = std::acos(cosTheta);
          y1 = arma::dot(a1, dR);
          y2 = arma::dot(a2, dR);
          sin2Theta = 1 - std::pow(cosTheta, 2);
          axis_shift_1 = (y1 + y2 * cosTheta) / sin2Theta;
          axis_shift_2 = (y2 + y1 * cosTheta) / sin2Theta;
          z_shift = arma::norm((axis_shift_1 * a1 + s1.pos()) -
                                (axis_shift_2 * a2 + s2.pos()));
        }

        double rate =
          scat_tab->get_rate(theta, z_shift, axis_shift_1, axis_shift_2);

        std::pair<double, scatterer*> p ={rate, &s2};

        neighbors_list.push_back(p);
      }
    };

    for (auto& l : close_scats) {
      for (auto& s: (*l)){
        check_neighbor(*this, *s);
      }
    }

    for (unsigned i=1; i<neighbors_list.size(); ++i){
      neighbors_list[i].first += neighbors_list[i-1].first;
    }

    _max_rate = neighbors_list.back().first;
    _inverse_max_rate = 1/_max_rate;

    return neighbors_list;
  };
} // mc namespace
