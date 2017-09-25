#include <iostream>
#include <memory>
#include <experimental/filesystem>
#include <ctime>
#include <cmath>
#include <fstream>

#include "monte_carlo.h"
// #include "particle.h"
#include "utility.h"
#include "gas_ff.h"
#include "forster_ff.h"

namespace mc
{
  monte_carlo::monte_carlo() // constructor
	{
		mc::init_random_number_generator();

		_temperature = 300;
		_beta = 1./(mc::kB*_temperature);
		_time = 0.;

		// set simulation geometry parameters, the rest are set automatically
		mc::t_float contact_length = 100.e-9;
		mc::t_float bulk_length = 1000.e-9;
		mc::t_float cross_section_1 = 100.e-9;
		mc::t_float cross_section_2 = 100.e-9;
		mc::t_float total_length = 2.*contact_length+bulk_length;

		// set simulation domain for boundary reflection
		_domain.first = {0., 0., 0.};
		_domain.second = {cross_section_1, cross_section_2, total_length};

		// upper and lower corners of the contact and bulk
		mc::arr1d contact_1_lower_corner = {0., 0., 0.};
		mc::arr1d contact_1_upper_corner = {cross_section_1, cross_section_2, contact_length};
		mc::arr1d bulk_lower_corner = {0., 0., contact_length};
		mc::arr1d bulk_upper_corner = {cross_section_1, cross_section_2, contact_length+bulk_length};
		mc::arr1d contact_2_lower_corner = {0., 0., contact_length+bulk_length};
		mc::arr1d contact_2_upper_corner = {cross_section_1, cross_section_2, total_length};

		_contacts.emplace_back(contact_1_lower_corner, contact_1_upper_corner);
		_contacts.emplace_back(contact_2_lower_corner, contact_2_upper_corner);
		_bulk.set_borders(bulk_lower_corner, bulk_upper_corner);

    _pilot = std::make_shared<mc::gas_free_flight>();
		_scatterer = std::make_shared<mc::scatter>();

		_contacts[0].populate(_beta, 1100, _pilot, _scatterer);
		_contacts[1].populate(_beta, 100, _pilot, _scatterer);
		_bulk.populate(_beta,0, _pilot, _scatterer);

		_number_of_profile_sections = 10;
		_population_profile.first = 0;
		_population_profile.second = std::vector<mc::t_uint>(_number_of_profile_sections, 0);
		_current_profile.first = 0;
		_current_profile.second = std::vector<mc::t_float>(_number_of_profile_sections, 0.);
	};
} // namespace mc
