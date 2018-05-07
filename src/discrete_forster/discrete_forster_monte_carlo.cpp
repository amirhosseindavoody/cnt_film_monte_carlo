#include <iostream>
#include <list>
#include <experimental/filesystem>
#include <fstream>
#include <map>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <armadillo>


#include "../helper/prepare_directory.hpp"
#include "../helper/progress.hpp"
#include "../exciton_transfer/cnt.h"
#include "discrete_forster_monte_carlo.h"


namespace mc
{

  // high level method to calculate proper scattering table
  void discrete_forster_monte_carlo::initialize_scattering_table() {
    std::cout << "\ninitializing scattering table...\n";

    switch(_scat_t) {
      case davoody:
        {
          // get the parent directory for cnts
          std::string parent_directory = _json_scat["cnts"]["directory"];
          _json_scat["cnts"].erase("directory");

          // create excitons and calculate exciton dispersions
          std::vector<cnt> cnts;
          cnts.reserve(_json_scat["cnts"].size()); // this is reservation of space is crucial to ensure we do not move cnts, since the move constructor is not implemented yet
          for (const auto& j_cnt: _json_scat["cnts"])
          {
            cnts.emplace_back(cnt(j_cnt,parent_directory));
            cnts.back().calculate_exciton_dispersion();
          };
          _scat_table = create_davoody_scatt_table(cnts[0], cnts[0]);
        }
        break;

      case forster:
        _scat_table = create_forster_scatt_table(1.e15, 1.4e9);
        break;

      case wong:
        _scat_table = create_forster_scatt_table(1.e13, 1.4e9);
        break;

      default:
        throw std::logic_error("invalid value for rate type");
    }
  };


  // method to calculate scattering rate via davoody et al. method
  discrete_forster_monte_carlo::scattering_struct discrete_forster_monte_carlo::create_davoody_scatt_table(const cnt& d_cnt, const cnt& a_cnt) {
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
  discrete_forster_monte_carlo::scattering_struct discrete_forster_monte_carlo::create_forster_scatt_table(double gamma_0, double r_0) {
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

  // find the neighbors of each scattering object
	void discrete_forster_monte_carlo::find_neighbors(std::list<std::shared_ptr<discrete_forster_scatter>>& scat_list, const double& max_hopping_radius, const double& max_search_radius) {

    std::cout << "finding neighbors in scatterers list:\n";

		// sort scattering objects to have better performance
		scat_list.sort([](const auto& s1, const auto& s2){
			return s1->pos(1) < s2->pos(1);
		});


		int counter = 0;
		double avg_max_rate = 0;
		double avg_number_of_neighbors = 0;
		double another_counter = 0;
    long total_number_of_neighbors=0;

		for (auto i = scat_list.begin(); i != scat_list.end(); ++i)
		{

			for (auto j = std::next(i); j!= scat_list.end(); ++j)
			{
        arma::vec dR = (*i)->pos()-(*j)->pos();
        double distance = arma::norm(dR);

				if ((distance < max_hopping_radius) and (distance > 0.4e-9))
				{
          arma::vec a1 = (*i)->orientation();
          arma::vec a2 = (*j)->orientation();
          double cosTheta = arma::dot(a1,a2);
          double theta = std::acos(cosTheta);
          double y1 = arma::dot(a1,dR);
          double y2 = arma::dot(a2,dR);
          double sin2Theta = 1-std::pow(cosTheta,2);
          double axis_shift_1 = (y1+y2*cosTheta)/sin2Theta;
          double axis_shift_2 = (y2+y1*cosTheta)/sin2Theta;
          double z_shift = arma::norm((axis_shift_1*a1+(*i)->pos())-(axis_shift_2*a2+(*j)->pos()));
          // check if parallel case has happend
          if (cosTheta == 1){
            axis_shift_1=0;
            axis_shift_2 = arma::dot(dR,a1);
            theta = 0;
            z_shift = arma::norm(dR-arma::dot(dR,a1)*a1);
          }

          double rate = _scat_table.get_rate(theta,z_shift,axis_shift_1,axis_shift_2);

					(*i)->add_neighbor(*j,distance,rate);
					(*j)->add_neighbor(*i,distance,rate);
          total_number_of_neighbors++;
				}

				// break the search if the scatterers are getting too far apart to speed up the search process
				if (distance > max_search_radius)
				{
					break;
				}
			}
			(*i)->sort_neighbors();
			(*i)->make_cumulative_scat_rate();

			counter ++;
			avg_number_of_neighbors += (*i)->number_of_neighbors();
			avg_max_rate += (*i)->max_rate();
			another_counter += 1.;

			if (counter %1000 == 0)
			{
				std::cout << "scatterer number:" << counter << " ... average number of neighbors:" << (avg_number_of_neighbors/another_counter) << " ... average max rate = " << avg_max_rate/another_counter << "                     \r" << std::flush;
				avg_max_rate = 0;
				avg_number_of_neighbors = 0;
				another_counter = 0;
			}

		}

		std::cout << "\n\nfinished finding neighbors:\n" 
              << "total number of neighbor pairs: " << total_number_of_neighbors << "\n"
              << "average number of neighbors per scatterer: " << (2*total_number_of_neighbors)/scat_list.size() << "\n\n";
	};

  // create a crystalline mesh structure
	std::list<std::shared_ptr<discrete_forster_scatter>> discrete_forster_monte_carlo::create_crystalline_structure()
  {
    std::cout << "\n\nmaking crystalline mesh of scatterer objects...\n";

    std::list<std::shared_ptr<discrete_forster_scatter>> scatterer_list;

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
        n_tubes(i) = unsigned(domain_length(i)/cnt_length);
        dr(i) = cnt_length;
      } else {
        n_tubes(i) = unsigned(domain_length(i)/cnt_spacing);
        dr(i) = cnt_spacing;
      }
    }
    n_tubes.print("number of tubes along each dimension:");

    for (unsigned ix=0; ix<=n_tubes(0); ix++){
      double x = ix*dr(0);
      for (unsigned iy=0; iy<=n_tubes(1); iy++){
        double y = iy*dr(1);
        for (unsigned iz=0; iz<=n_tubes(2); iz++){
          double z = iz*dr(2);

          scatterer_list.push_back(std::make_shared<discrete_forster_scatter>());
          scatterer_list.back()->set_pos({x,y,z});
          scatterer_list.back()->set_orientation(orientation);
        }
      }
    }

    return scatterer_list;
  };

  // find the neighbors of each scattering object using segmentation method,
  // first divide the scattering objects into multiple buckets based on their
  // position in x-z plane, and then search for the neighbor pairs by search only
  // through the neighboring buckets.
  void discrete_forster_monte_carlo::find_neighbors_using_bucket(
      std::list<std::shared_ptr<discrete_forster_scatter>>& scat_list,
      const double max_hopping_radius) {
    std::cout << "finding neighbors in scatterers list using buckets:\n";

    double xmin = (_domain.first)(0);
    double xmax = (_domain.second)(0);
    int nx = std::ceil((xmax - xmin)/max_hopping_radius);

    double zmin = (_domain.first)(2);
    double zmax = (_domain.second)(2);
    int nz = std::ceil((zmax - zmin) / max_hopping_radius)+1;


    std::cout << "xmin:" << xmin << " , xmax:" << xmax << "\n";
    std::cout << "zmin:" << zmin << " , zmax:" << zmax << "\n";
    std::cout << "nx:" << nx << " , nz:" << nz << "\n";

    typedef std::shared_ptr<discrete_forster_scatter> s_ptr;
    std::vector<std::vector<std::list<s_ptr>>> grid(nx, std::vector<std::list<s_ptr>>(nz, std::list<s_ptr>()));

    for (s_ptr& s: scat_list){
      int ix = (s->pos(0)-xmin)/max_hopping_radius;
      int iz = (s->pos(2)-zmin)/max_hopping_radius;
      grid[ix][iz].push_back(s);
    }

    // a local reference to the scattering table since we cannot access class members via capture by reference in lambda function
    scattering_struct& scat_tab_ptr = _scat_table;
    auto find_pairs = [&max_hopping_radius, &scat_tab_ptr](std::list<s_ptr>& l1, std::list<s_ptr>& l2) {
      double cosTheta, theta, y1, y2, sin2Theta, axis_shift_1, axis_shift_2, z_shift;

      for (s_ptr& s1: l1){
        for (s_ptr& s2: l2){
          arma::vec dR = s1->pos() - s2->pos();
          double distance = arma::norm(dR);

          if ((distance < max_hopping_radius) and (distance > 0.4e-9)) {
            arma::vec a1 = s1->orientation();
            arma::vec a2 = s2->orientation();
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
              z_shift = arma::norm((axis_shift_1 * a1 + s1->pos()) -
                                          (axis_shift_2 * a2 + s2->pos()));
            }

            double rate = scat_tab_ptr.get_rate(theta, z_shift, axis_shift_1,
                                               axis_shift_2);

            s1->add_neighbor(s2, distance, rate);
          }

        }
      }
    };



    for (int ix=0; ix<nx; ++ix){
      std::cout << "ix:" << ix << "\r" << std::flush;
      for (int iz=0; iz<nz; ++iz){
        find_pairs(grid[ix][iz], grid[ix][iz]);
        
        if (ix > 0) find_pairs(grid[ix][iz], grid[ix - 1][iz]);
        if (ix + 1 < nx) find_pairs(grid[ix][iz], grid[ix + 1][iz]);
        if (iz > 0) find_pairs(grid[ix][iz], grid[ix][iz - 1]);
        if (iz + 1 < nz) find_pairs(grid[ix][iz], grid[ix][iz + 1]);

        if (ix > 0 && iz > 0) find_pairs(grid[ix][iz], grid[ix - 1][iz - 1]);
        if (ix + 1 < nx && iz > 0) find_pairs(grid[ix][iz], grid[ix + 1][iz - 1]);
        if (ix > 0 && iz + 1 < nz) find_pairs(grid[ix][iz], grid[ix - 1][iz + 1]);
        if (ix + 1 < nx && iz + 1 < nz) find_pairs(grid[ix][iz], grid[ix + 1][iz + 1]);
      }
    }

    
    
    int counter = 0;
    double avg_max_rate = 0;
    double avg_number_of_neighbors = 0;
    double another_counter = 0;
    long total_number_of_neighbors = 0;

    for (s_ptr& s:scat_list) {
      
      s->sort_neighbors();
      s->make_cumulative_scat_rate();

      counter++;
      avg_number_of_neighbors += s->number_of_neighbors();
      avg_max_rate += s->max_rate();
      another_counter += 1.;

      if (counter % 1000 == 0) {
        std::cout << "scatterer number:" << counter
                  << " ... average number of neighbors:"
                  << (avg_number_of_neighbors / another_counter)
                  << " ... average max rate = "
                  << avg_max_rate / another_counter << "                     \r"
                  << std::flush;
        avg_max_rate = 0;
        avg_number_of_neighbors = 0;
        another_counter = 0;
      }
    }

    std::cout << "\n\nfinished finding neighbors:\n"
              << "total number of neighbor pairs: " << total_number_of_neighbors
              << "\n"
              << "average number of neighbors per scatterer: "
              << (2 * total_number_of_neighbors) / scat_list.size() << "\n\n";
  };

} // end of namespace mc