#include <iostream>
#include <list>
#include <experimental/filesystem>
#include <fstream>
#include <map>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <armadillo>

#include "../../lib/json.hpp"
#include "../exciton_transfer/cnt.h"
#include "../helper/prepare_directory.hpp"
#include "../helper/progress.hpp"
#include "monte_carlo.h"

namespace mc
{

  // high level method to calculate proper scattering table
  scattering_struct monte_carlo::create_scattering_table(nlohmann::json j) {
    if (j.count("rate type")==0)
      throw std::invalid_argument(
          "\"rate type\" must be specifieced!!");

    std::cout << "\ninitializing scattering table...\n";

    scattering_struct scat_tab;

    if (j["rate type"].get<std::string>() == "davoody") {

      // get the parent directory for cnts
      std::string parent_directory = j["cnts"]["directory"];
      j["cnts"].erase("directory");

      // create excitons and calculate exciton dispersions
      std::vector<cnt> cnts;
      cnts.reserve(j["cnts"].size()); // this is reservation of space is crucial to ensure we do not move
                                      // cnts, since the move constructor is not implemented yet
      
      for (const auto& j_cnt : j["cnts"]) {
        cnts.emplace_back(cnt(j_cnt, parent_directory));
        cnts.back().calculate_exciton_dispersion();
      };
      scat_tab = create_davoody_scatt_table(cnts[0], cnts[0]);


    } else if (j["rate type"].get<std::string>() == "forster") {
      scat_tab = create_forster_scatt_table(1.e15, 1.4e9);
    } else if (j["rate type"].get<std::string>() == "wong") {
      scat_tab = create_forster_scatt_table(1.e13, 1.4e9);
    } else {
      throw std::invalid_argument(
          "rate type must be one of the following: \"davoody\", \"forster\", \"wong\"");
    }

    return scat_tab;
  };

  // method to calculate scattering rate via davoody et al. method
  scattering_struct monte_carlo::create_davoody_scatt_table(const cnt& d_cnt, const cnt& a_cnt) {
    auto zshift_prop = _json_prop["zshift [m]"];
    arma::vec z_shift = arma::linspace<arma::vec>(zshift_prop[0], zshift_prop[1], zshift_prop[2]);

    auto axis_shift_prop_1 = _json_prop["axis shift 1 [m]"];
    arma::vec axis_shift_1 = arma::linspace<arma::vec>(axis_shift_prop_1[0], axis_shift_prop_1[1], axis_shift_prop_1[2]);

    auto axis_shift_prop_2 = _json_prop["axis shift 2 [m]"];
    arma::vec axis_shift_2 = arma::linspace<arma::vec>(axis_shift_prop_2[0], axis_shift_prop_2[1], axis_shift_prop_2[2]);

    auto theta_prop = _json_prop["theta [degrees]"];
    arma::vec theta = arma::linspace<arma::vec>(theta_prop[0], theta_prop[1], theta_prop[2])*(constants::pi/180);

    arma::field<arma::cube> rate(theta.n_elem);
    rate.for_each([&](arma::cube& c){c.zeros(z_shift.n_elem, axis_shift_1.n_elem, axis_shift_2.n_elem);});

    exciton_transfer ex_transfer(d_cnt, a_cnt);

    progress_bar prog(theta.n_elem*z_shift.n_elem*axis_shift_1.n_elem*axis_shift_2.n_elem,"create davoody scattering table");

    double max_rate = 0;
    double min_rate = 10e15;

    unsigned i_th=0;
    for (const auto& th: theta)
    {
      unsigned i_zsh=0;
      for (const auto& zsh: z_shift)
      {
        unsigned i_ash1=0;
        for (const auto& ash1: axis_shift_1)
        {
          unsigned i_ash2=0;
          for (const auto& ash2: axis_shift_2)
          {
            prog.step();
            rate(i_th)(i_zsh,i_ash1,i_ash2) = ex_transfer.first_order(zsh, {ash1, ash2}, th, false);

            if (rate(i_th)(i_zsh,i_ash1,i_ash2) > max_rate) {
              max_rate = rate(i_th)(i_zsh,i_ash1,i_ash2);
            }
            if (rate(i_th)(i_zsh,i_ash1,i_ash2) < min_rate) {
              min_rate = rate(i_th)(i_zsh,i_ash1,i_ash2);
            }

            i_ash2++;
          }
          i_ash1++;
        }
        i_zsh++;
      }
      i_th++;
    }

    scattering_struct scat_table(rate,theta,z_shift,axis_shift_1,axis_shift_2);

    std::cout << "\nmax rate in davoody scattering table: " << max_rate << " [1/s]\n";
    std::cout << "min rate in davoody scattering table: " << min_rate << " [1/s]\n\n";

    return scat_table;
  };

  // method to calculate scattering rate via forster method
  scattering_struct monte_carlo::create_forster_scatt_table(double gamma_0, double r_0) {
    auto zshift_prop = _json_prop["zshift [m]"];
    arma::vec z_shift = arma::linspace<arma::vec>(zshift_prop[0], zshift_prop[1], zshift_prop[2]);

    auto axis_shift_prop_1 = _json_prop["axis shift 1 [m]"];
    arma::vec axis_shift_1 = arma::linspace<arma::vec>(axis_shift_prop_1[0], axis_shift_prop_1[1], axis_shift_prop_1[2]);

    auto axis_shift_prop_2 = _json_prop["axis shift 2 [m]"];
    arma::vec axis_shift_2 = arma::linspace<arma::vec>(axis_shift_prop_2[0], axis_shift_prop_2[1], axis_shift_prop_2[2]);

    auto theta_prop = _json_prop["theta [degrees]"];
    arma::vec theta = arma::linspace<arma::vec>(theta_prop[0], theta_prop[1], theta_prop[2])*(constants::pi/180);

    arma::field<arma::cube> rate(theta.n_elem);
    rate.for_each([&](arma::cube& c){c.zeros(z_shift.n_elem, axis_shift_1.n_elem, axis_shift_2.n_elem);});

    progress_bar prog(theta.n_elem*z_shift.n_elem*axis_shift_1.n_elem*axis_shift_2.n_elem,"create forster scattering table");

    unsigned i_th=0;
    for (const auto& th: theta)
    {
      unsigned i_zsh=0;
      for (const auto& zsh: z_shift)
      {
        unsigned i_ash1=0;
        for (const auto& ash1: axis_shift_1)
        {
          unsigned i_ash2=0;
          for (const auto& ash2: axis_shift_2)
          {
            prog.step();
            arma::vec r1 = {ash1, 0, 0};
            arma::vec r2 = {ash2*std::cos(th), ash2*std::cos(th), zsh};
            arma::vec dR = r1-r2;
            double angle_factor = std::cos(th)-3*arma::dot(arma::normalise(r1),arma::normalise(dR))*arma::dot(arma::normalise(r2),arma::normalise(dR));
            rate(i_th)(i_zsh,i_ash1,i_ash2) = gamma_0*std::pow(angle_factor,2)*std::pow(1.e-9/arma::norm(dR),6);
            i_ash2++;
          }
          i_ash1++;
        }
        i_zsh++;
      }
      i_th++;
    }

    scattering_struct scat_table(rate,theta,z_shift,axis_shift_1,axis_shift_2);

    return scat_table;
  };

  // create a crystalline mesh structure
	std::vector<scatterer> monte_carlo::create_crystalline_structure() {
    std::cout << "\n\nmaking crystalline mesh of scatterer objects...\n";

    std::vector<scatterer> scatterer_list;


    arma::vec domain_length;
    arma::vec orientation;
    double cnt_length;
    double cnt_spacing;

    if (_json_prop.count("crystall mesh properties")){
      cnt_length = _json_prop["crystall mesh properties"]["cnt length [m]"];
      cnt_spacing = _json_prop["crystall mesh properties"]["cnt spacing [m]"];

      std::vector<double> tmp1 = _json_prop["crystall mesh properties"]["cnt orientation"];
      orientation = tmp1;

      std::vector<double> tmp2 = _json_prop["crystall mesh properties"]["simulation domain dimensions [m]"];
      domain_length = tmp2;
    } else {
      throw std::invalid_argument("must specify \"crystall mesh properties\" in input.json");
    }

    orientation = arma::normalise(orientation);

    int dim = 0;
    for (;dim<3;dim++){
      if (orientation(dim)!=0){
        break;
      }
    }
    std::cout << "orientation of the tubes are along " << dim << "'th axis\n";

    arma::uvec n_tubes = {0,0,0};
    arma::vec dr = {0,0,0};
    for (int i=0; i<3; i++){
      if (i==dim){
        n_tubes(i) = unsigned(domain_length(i)/cnt_length)+1;
        dr(i) = cnt_length;
      } else {
        n_tubes(i) = unsigned(domain_length(i)/cnt_spacing)+1;
        dr(i) = cnt_spacing;
      }
    }
    n_tubes.print("number of tubes along each dimension:");

    scatterer_list.resize(n_tubes(0)*n_tubes(1)*n_tubes(2));

    for (unsigned ix=0; ix<n_tubes(0); ix++){
      double x = ix*dr(0);
      for (unsigned iy=0; iy<n_tubes(1); iy++){
        double y = iy*dr(1);
        for (unsigned iz=0; iz<n_tubes(2); iz++){
          double z = iz*dr(2);

          unsigned idx = ix + iy*n_tubes(0) + iz * n_tubes(0) * n_tubes(1);
          scatterer_list[idx].set_pos({x,y,z});
          scatterer_list[idx].set_orientation({orientation(0), orientation(1), orientation(2)});
        }
      }
    }

    return scatterer_list;
  };

} // end of namespace mc