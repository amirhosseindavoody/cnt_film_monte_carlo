#include <iostream>
#include <memory>
#include <experimental/filesystem>
#include <ctime>
#include <cmath>
#include <fstream>
#include <cstring>

#include "../rapidxml/rapidxml.hpp"
#include "../rapidxml/rapidxml_utils.hpp"
#include "../rapidxml/rapidxml_print.hpp"
#include "monte_carlo.h"
#include "utility.h"


namespace mc
{
  // constructor
  monte_carlo::monte_carlo()
	{
		mc::init_random_number_generator();

		set_temperature(300);
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

    _pilot = std::make_shared<mc::t_free_flight>();
    // _scatterer = std::make_shared<mc::t_scatter>();
    _scatterer = std::make_shared<mc::t_scatter>();

		_contacts[0].populate(_beta, 1100, _pilot, _scatterer);
		_contacts[1].populate(_beta, 100, _pilot, _scatterer);
		_bulk.populate(_beta,0, _pilot, _scatterer);

		_number_of_profile_sections = 10;
		_population_profile.first = 0;
		_population_profile.second = std::vector<mc::t_uint>(_number_of_profile_sections, 0);
		_current_profile.first = 0;
		_current_profile.second = std::vector<mc::t_float>(_number_of_profile_sections, 0.);
	};

  // set the output directory and the output file name
  void monte_carlo::process_command_line_args(int argc, char* argv[])
  {
    namespace fs = std::experimental::filesystem;

    fs::directory_entry xml_file;

    std::cout << "current path is " << fs::current_path() << std::endl;

    if (argc <= 1)
    {
      xml_file.assign("input.xml");
    }
    else
    {
      xml_file.assign(argv[1]);
    }

    if(fs::exists(xml_file))
    {
      std::cout << "input xml file found: " << xml_file.path() << std::endl;
    }
    else
    {
      std::cout << "input xml file NOT found: " << xml_file.path() << std::endl;
      std::exit(1);
    }

    if (!fs::is_regular_file(xml_file))
    {
      std::cout << "input xml file NOT found: " << xml_file.path() << std::endl;
      std::exit(1);
    }
    std::cout << std::endl;

    rapidxml::file<> xmlFile(xml_file.path().c_str()); //open file
    rapidxml::xml_document<> doc; //create xml object
    doc.parse<0>(xmlFile.data()); //parse contents of file
    rapidxml::xml_node<>* curr_node = doc.first_node(); //gets the node "Document" or the root nodes

    // set the output_directory
    {
      curr_node = curr_node->first_node("output_directory");
      std::string attr = curr_node->first_attribute("type")->value();
      std::string path = mc::trim(curr_node->value());
      if (attr == "absolute")
      {
        std::cout << "absolute directory format used!\n";
      }

      _output_directory.assign(path);
      std::cout << "output_directory: " << _output_directory.path() << std::endl;

      if (not fs::exists(_output_directory.path()))
      {
        std::cout << "warning: output directory does NOT exist!!!" << std::endl;
        std::cout << "output directory: " << _output_directory.path() << std::endl;
        fs::create_directories(_output_directory.path());
      }

      if (fs::is_directory(_output_directory.path()))
      {
        if (not fs::is_empty(_output_directory.path()))
        {
          std::cout << "warning: output directory is NOT empty!!!" << std::endl;
          std::cout << "output directory: " << _output_directory.path() << std::endl;
          std::cout << "deleting the existing directory!!!" << std::endl;
          fs::remove_all(_output_directory.path());
          fs::create_directories(_output_directory.path());
        }
      }
      else
      {
        std::cout << "error: output path is NOT a directory!!!" << std::endl;
        std::cout << "output path: " << _output_directory.path() << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // set the temperature
    {
      auto next_node = curr_node->next_sibling("temperature");
      if (next_node == 0)
      {
        next_node = curr_node->previous_sibling("temperature");
        if (next_node == 0)
        {
          std::cout << "temperature not found in XML file!!!" << std::endl;
          std::exit(1);
        }
      }
      curr_node = next_node;

      std::string attr = curr_node->first_attribute("units")->value();

      if (attr == "Kelvin" or attr == "kelvin" or attr == "K" or attr == "k")
      {
        set_temperature(std::atof(curr_node->value()));
      }
      else if (attr == "celcius" or attr == "Celcius" or attr == "c" or attr == "C")
      {
        set_temperature(273.15 + std::atof(curr_node->value()));
      }
      else
      {
        std::cout << "error in temperature units!!!" << std::endl;
        std::exit(1);
      }
    }

    // set input directory
    {
      auto next_node = curr_node->next_sibling("input_directory");
      if (next_node == 0)
      {
        next_node = curr_node->previous_sibling("input_directory");
        if (next_node == 0)
        {
          std::cout << "input_directory not found in XML file!!!" << std::endl;
          std::exit(1);
        }
      }
      curr_node = next_node;

      std::string attr = curr_node->first_attribute("type")->value();
      std::string path = mc::trim(curr_node->value());
      if (attr == "absolute")
      {
        std::cout << "absolute directory format used!\n";
      }
      _input_directory.assign(path);

      if (not fs::exists(_input_directory.path()))
      {
        std::cout << "\n***\nwarning: input directory does NOT exist!!!\n"
                  << "input directory: " << _input_directory.path() << "\n***\n\n";
        std::exit(1);
      }

      if (fs::is_directory(_input_directory.path()))
      {
        if (fs::is_empty(_input_directory.path()))
        {
          std::cout << "warning: input directory is empty!!!" << std::endl;
          std::cout << "input directory: " << _input_directory.path() << std::endl;
        }
      }
      else
      {
        std::cout << "\n***\nerror: input path is NOT a directory!!!\n"
                  << "input path: " << _input_directory.path() << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // std::exit(0);

  };


} // namespace mc
