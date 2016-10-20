/**
CNT.cpp
Stores all relevant information for a carbon nanotube

@author Alex Gabourie
@version 1.00
*/

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <regex>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "CNT.h"
#include "write_log.h"

//Global variable
extern double ymax;

/**
Sets the CNT object to some default values. DO NOT USE CNTs CONSTRUCTED
THIS WAY. This is only to appease the compiling gods.

@return CNT Object
*/
CNT::CNT()
{
	n = 0;
	m = 0;
	length = 0;
	cylinderHeight = 0;
	tubeSeparation = 0;
	minSpacing = 0;
	diameter = 0;
	cntNum = 0;
	number_of_points = 0;
	positions = vector<vector<double>>(3);

}
/**
Reads a CNT file and creates a CNT object with all the information stored
in that file.

@param filePath The path of the file containing the CNT info
@return CNT Object
*/
CNT::CNT(const string fileName, const string folderPath, double segment_length)
{
	string filePath = folderPath + "/" + fileName;

	regex rgx("\\d+"); //basic_regex instantiation of type char
	{
		//Extract the tube number from the file path
		//filePath is the target sequence
		smatch matches; //match_results for string objects
		//search to see if sequence matches any part of target sequence
		regex_search(fileName, matches, rgx);
		//input checking
		if (!matches.empty())
		{
			cntNum = stoi(matches[0]);
		}
		//Cannot extract CNT number from file name.
		else
		{
			cout << "Error: Incorrect file names. File names should look like \"CNT_Num_x.csv\"" << " where\n \'x\' is a decimal number.\n";
			exit(EXIT_FAILURE);
		}
	}
	

	//Open .csv file to read.
	ifstream file(filePath);
	//Checks if file is open and readable
	if (!file.good())
	{
		cout << "Cannot read " << fileName << endl;
		exit(EXIT_FAILURE);
	}
	//Can now read the file
	string temp = " "; //stores intermediate strings during parsing

	//Chirality, first line//
	getline(file, temp, '\n');
	{
		std::string::size_type n = temp.find(",");
		if (n == std::string::npos)
		{
			cout << "Error in reading CNT chirality!!!" << endl;
			exit(EXIT_FAILURE);
		}

		temp = temp.substr(n+1);
		temp.erase(temp.size()-1);
		rgx.assign("(\\d+)(\\s+)(\\d+)");
		smatch matches; //match_results for string objects
		regex_match(temp, matches, rgx);
		if (!matches.empty())
		{
			n = stoi(matches[1]);
			m = stoi(matches[3]);
		}
		else
		{
			cout << "Error: Cannot extract chirality!!!" << endl;
			exit(EXIT_FAILURE);
		}
	}

	//CNT Length, second line
	getline(file, temp, '\n');
	{
		std::string::size_type n = temp.find(",");
		if (n == std::string::npos)
		{
			cout << "Error in reading CNT chirality!!!" << endl;
			exit(EXIT_FAILURE);
		}
		temp = temp.substr(n+1);
		
		length = stod(temp);
	}

	//Cylinder height, third line
	getline(file, temp, '\n');
	{
		std::string::size_type n = temp.find(",");
		if (n == std::string::npos)
		{
			cout << "Error in reading CNT chirality!!!" << endl;
			exit(EXIT_FAILURE);
		}
		temp = temp.substr(n+1);

		cylinderHeight = stod(temp);
	}

	//Intertube spacing, fourth line
	getline(file, temp, '\n');
	{
		std::string::size_type n = temp.find(",");
		if (n == std::string::npos)
		{
			cout << "Error in reading CNT chirality!!!" << endl;
			exit(EXIT_FAILURE);
		}
		temp = temp.substr(n+1);

		minSpacing = stod(temp);

	}
	
	//intercylinder spacing, fifth line
	getline(file, temp, '\n');
	{
		std::string::size_type n = temp.find(",");
		if (n == std::string::npos)
		{
			cout << "Error in reading CNT chirality!!!" << endl;
			exit(EXIT_FAILURE);
		}
		temp = temp.substr(n+1);

		tubeSeparation = stod(temp);

	}
	
	//Now that all parameters are extracted, calculate diameter
	setDiameter(n, m);

	//move past coordinate labels line
	getline(file, temp, '\n');

	//Counting line loop
	number_of_points = -1; //count of positions, negative b/c while loop exits
						//after going an extra line
	while (!file.eof())
	{
		getline(file, temp,'\n');
		number_of_points++;
	} 

	positions = vector<vector<double>>(3);
	for (int i = 0; i < 3; i++)
	{
		positions[i].resize(number_of_points);
	}

	file.clear(); //reset state flags
	file.seekg(0); //move to beginning of file
	//move to correct line to start position parsing
	for (int i = 0; i < 6; i++)
	{
		getline(file, temp, '\n');
	}

	
	for (int i = 0; i < number_of_points && file.good(); i++)
	{

		string tmp_string = " ";

		getline(file, temp, '\n');
		
		if (file.eof())
		{
			break;
		}
		
		//Grab all of the position data
		istringstream ss(temp);
		getline(ss, tmp_string, ',');
		double x = stod(tmp_string);
		
		getline(ss, tmp_string, ',');
		double y = stod(tmp_string);
		if (y > ymax){ ymax = y; }
		
		getline(ss, tmp_string, ',');
		double z = stod(tmp_string);

		positions[0][i] = x;
		positions[1][i] = y;
		positions[2][i] = z;
	}
	
	//Calculate the segments needed for table generation
	calculate_segments(segment_length);

	//Checks to see if some segments were calculated
	if (segments->empty())
	{
		cout << "Error: No segments calculated for tube number: ";
		cout << cntNum;
		cout << "\n";
		exit(EXIT_FAILURE);
	}

	initialized = true;

}

//Destructor implemented automatically


/**
Uses the chirality of the nanotube to calculate the diameter
@param n Hamada n paramter
@param m Hamada m parameter
*/
void CNT::setDiameter(int n, int m)
{
	diameter = A_CC*sqrt(pow(n*1., 2) + pow(m*1., 2) + n*m) / M_PI;
}

/**
Gets the diameter of the CNT

@param void
@return Diameter of the CNT
*/
double CNT::getDiameter()
{
	return diameter;
}

/**
Gets the length of the CNT

@param void
@return Length of the CNT
*/
double CNT::getLength()
{
	return length;
}

/**
Gets the height of each compositional cylinders

@param void
@return The height of each compositional cylinders
*/
double CNT::getCylinderHeight()
{
	return cylinderHeight;
}

/**
Gets the separation between two compositional cylinders

@param void
@return The separation between two compositional cylinders
*/
double CNT::getTubeSeparation()
{
	return tubeSeparation;
}

/**
Gets the minimum spacing between two nanotubes

@param void
@return minimum spacing between two nanotubes
*/
double CNT::getMinSpacing()
{
	return minSpacing;
}

/**
Gets m parameter of the CNT

@param void
@return m
*/
int CNT::getm()
{
	return m;
}

/**
Gets n parameter of the CNT

@param void
@return n
*/
int CNT::getn()
{
	return n;
}


/**
Gets the tube number

@param void
@return CNT number
*/
int CNT::getCNTNum()
{
	return cntNum;
}

/**
Says whether or not the CNT was initialized

@return initialization status
*/
bool CNT::isInitialized()
{
	return initialized;
}
 
/**
Calculates the segments used for the MC simulations

@param segment_length The desired length of the segments
@return Vector of the segments
*/
void CNT::calculate_segments(double segment_length)
{
	//parameter check
	if (length < segment_length)
	{
		cout << "Error: Tube length is smaller than minimum segment length!!!" << endl;
		exit(EXIT_FAILURE);
	}
	else if (segment_length <= 0)
	{
		cout << "Error: Minimum segment length must be positive!!!" << endl;
		exit(EXIT_FAILURE);
	}

	// segments = make_shared<vector<shared_ptr<segment>>>();

	double current_length = 0.;

	int first_idx=0;
	int second_idx=0;
	int segment_number = 0;

	for (int i = 0; i<number_of_points-1; i++)
	{

		double tmp_length = 0.0;
		for (int dim = 0; dim < 3; dim++)
		{
			double delta = positions[dim][i]-positions[dim][i+1];
			tmp_length += pow(delta,2);
		}
		tmp_length = sqrt(tmp_length);

		if (tmp_length > segment_length)
		{
			cout << "segment_length is too small for current CNT spine mesh!!!" << endl;
			exit(EXIT_FAILURE);
		}
		current_length += tmp_length;


		if (current_length > segment_length)
		{
			second_idx = i+1;
			int n = second_idx-first_idx;
			segment_number ++;

			vector<double> first_point(3);
			vector<double> second_point(3);

			perform_PCA(first_idx, n, segment_length, first_point, second_point);
			// print_segment_points(first_idx, n);

			segment curr_segment(segment_number, first_point, second_point);
			segments_new.push_back(curr_segment);

			first_idx = second_idx;
			current_length = 0.;
		}
	}

	exit(EXIT_SUCCESS);

}


void CNT::print_segment_points(int first_idx, int n)
{

	stringstream log_input;

	for (int j = 0; j<n; j++)
	{
		log_input 	<< std::scientific << positions[0][first_idx+j] << "    "
					<< std::scientific << positions[1][first_idx+j] << "    "
					<< std::scientific << positions[2][first_idx+j]
					<< endl;
	}

	ofstream my_file;
	my_file.open("bullet_physics_points.dat", ios::app);
	my_file << log_input.str();

	if (my_file.fail())
	{
		cout << "error in writing to a file!!!" << endl;
		exit(EXIT_FAILURE);
	}

	my_file.close();
}


void CNT::perform_PCA(int first_idx, int n, double segment_length, vector<double> &first_point, vector<double> &second_point)
{
	gsl_vector *avg = gsl_vector_alloc(3);
	gsl_vector_set_zero(avg);
	for (int j = 0; j<n; j++)
	{
		for (int dim=0; dim<3; dim++)
		{
			double val = gsl_vector_get(avg, dim) + positions[dim][first_idx+j];
			gsl_vector_set(avg, dim, val);
		}
	}
	gsl_vector_scale(avg, 1.0/double(n));
	
	gsl_matrix *X = gsl_matrix_alloc(n,3);
	for (int j=0; j<n; j++)
	{
		for (int dim=0; dim<3; dim++)
		{
			double val = positions[dim][first_idx+j] - gsl_vector_get(avg, dim);
			gsl_matrix_set(X, j, dim, val);

		}
	}

	gsl_matrix* XT = gsl_matrix_alloc(3, n);
	gsl_matrix_transpose_memcpy(XT, X);

	gsl_matrix* cov = gsl_matrix_alloc(3,3);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1/double(n), XT, X, 0.0, cov);


	gsl_vector *work = gsl_vector_alloc(3);
	gsl_vector *S = gsl_vector_alloc(3);
	gsl_matrix *V = gsl_matrix_alloc(3,3);
	gsl_linalg_SV_decomp (cov, V, S, work);

	cout << gsl_matrix_get(cov, 0, 0) << "   "
		 << gsl_matrix_get(cov, 1, 0) << "   "
		 << gsl_matrix_get(cov, 2, 0) << endl;

	cout << "norm = " << sqrt(pow(gsl_matrix_get(cov, 0, 0),2) + pow(gsl_matrix_get(cov, 1, 0),2) + pow(gsl_matrix_get(cov, 2, 0),2)) << endl;

	for (int dim=0; dim<3; dim++)
	{
		double val1 = +gsl_matrix_get(cov, dim, 0)*segment_length/2.0 + gsl_vector_get(avg, dim);
		double val2 = -gsl_matrix_get(cov, dim, 0)*segment_length/2.0 + gsl_vector_get(avg, dim);
		first_point[dim] = val1;
		second_point[dim] = val2;
	}

	// stringstream log_input;
	// log_input 	<< std::scientific << first_point[0] << "    "
	// 			<< std::scientific << first_point[1] << "    "
	// 			<< std::scientific << first_point[2]
	// 			<< endl
	// 		 	<< std::scientific << second_point[0] << "    "
	// 			<< std::scientific << second_point[1] << "    "
	// 			<< std::scientific << second_point[2]
	// 			<< endl;

	// ofstream my_file;
	// my_file.open("pca_points.dat", ios::app);
	// my_file << log_input.str();

	// if (my_file.fail())
	// {
	// 	cout << "error in writing to a file!!!" << endl;
	// 	exit(EXIT_FAILURE);
	// }

	// my_file.close();

	gsl_vector_free(work);
	gsl_vector_free(S);
	gsl_matrix_free(X);
	gsl_matrix_free(XT);
	gsl_matrix_free(cov);
	gsl_matrix_free(V);

}

