#include <iostream>
#include <limits>

#include "./particle.h"
#include "./scatterer.h"

namespace mc
{
  // update the final state of the particle
  const scatterer* scatterer::update_state(const double& max_hopping_radius) const {

    using namespace std;
    vector<pair<double, scatterer*>> nlist = find_neighbors(max_hopping_radius);

    if (nlist.empty())
      return this;

    double dice = nlist.back().first * double(std::rand()) / double(RAND_MAX);
    unsigned left = -1;
    unsigned right = nlist.size()-1;

    while (left + 1 < right) {
      unsigned mid = (left + right) / 2;
      if (nlist[mid].first <= dice) {
        left = mid;
      } else {
        right = mid;
      }
    }

    return nlist[right].second;
	};

  // find neighbors of the current scatterer and their scattering rates
  std::vector<std::pair<double, scatterer*>> scatterer::find_neighbors(const double& max_hopping_radius) const {
    std::vector<std::pair<double, scatterer*>> neighbors_list;

    auto check_neighbor = [&max_hopping_radius, &neighbors_list, this](scatterer& s2) {
      double cosTheta, theta, y1, y2, sin2Theta, axis_shift_1, axis_shift_2, z_shift;

      std::valarray<double> dR = this->pos() - s2.pos();
      double distance = std::sqrt((dR*dR).sum());

      if ((distance < max_hopping_radius) && (distance > 0.4e-9)) {
        std::valarray<double> a1 = this->orientation();
        std::valarray<double> a2 = s2.orientation();
        cosTheta = (a1*a2).sum();
        // check if parallel case has happend
        if (cosTheta == 1) {
          axis_shift_1 = 0;
          axis_shift_2 = (dR * a1).sum();
          theta = 0;
          z_shift = std::sqrt(std::pow(dR - (dR*a1).sum()*a1, 2).sum());
        } else {
          theta = std::acos(cosTheta);
          y1 = (a1*dR).sum();
          y2 = (a1*dR).sum();
          sin2Theta = 1 - std::pow(cosTheta, 2);
          axis_shift_1 = (y1 + y2 * cosTheta) / sin2Theta;
          axis_shift_2 = (y2 + y1 * cosTheta) / sin2Theta;
          z_shift = std::sqrt(std::pow((axis_shift_1 * a1 + this->pos()) - (axis_shift_2 * a2 + s2.pos()), 2).sum());
        }

        double rate = scat_tab->get_rate(theta, z_shift, axis_shift_1, axis_shift_2);

        std::pair<double, scatterer*> p = {rate, &s2};

        neighbors_list.push_back(p);
      }
    };

    for (auto& l : close_scats) {
      for (auto& s : (*l)) {
        check_neighbor(*s);
      }
    }

    // make accumulative scattering table
    for (unsigned i = 1; i < neighbors_list.size(); ++i) {
      neighbors_list[i].first += neighbors_list[i - 1].first;
    }

    return neighbors_list;
  };

  int scatterer::no_of_neighbors(const double& max_hopping_radius) const {
    int count=0;
    
    auto check_neighbor = [&max_hopping_radius, this, &count](scatterer* s2) {
      std::valarray<double> dR = this->pos() - s2->pos();
      double distance = std::sqrt((dR*dR).sum());

      if ((distance < max_hopping_radius) && (distance > 0.4e-9)) {
        count++;
      }
    };

    for (auto& l : close_scats) {
      for (auto& s : (*l)) {
        check_neighbor(s);
      }
    }

    return count;
  };
} // mc namespace
