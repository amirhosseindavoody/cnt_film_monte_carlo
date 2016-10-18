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

	coordinates = gsl_matrix_alloc(number_of_points, 3);

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

		gsl_matrix_set(coordinates, i, 0, x);
		gsl_matrix_set(coordinates, i, 1, y);
		gsl_matrix_set(coordinates, i, 2, z);

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
 
// /**
// Calculates the segments used for the MC simulations

// @param segment_length The desired length of the segments
// @return Vector of the segments
// */
// shared_ptr<vector<shared_ptr<segment>>> CNT::calculate_segments(double segLenMin)
// {
// 	//parameter check
// 	if (length < segLenMin)
// 	{
// 		cout << "Error: Tube length is smaller than minimum segment length!!!" << endl;
// 		exit(EXIT_FAILURE);
// 	}
// 	else if (segLenMin <= 0)
// 	{
// 		cout << "Error: Minimum segment length must be positive!!!" << endl;
// 		exit(EXIT_FAILURE);
// 	}

// 	//number of segments to use
// 	int numSegs = static_cast<int>( length / segLenMin); 
// 	//extra length past numSegs*segLenMin, important for expanding segment_length
// 	double extra = length - segLenMin*numSegs;
// 	//The equally lengthed segment lengths
// 	// half the length of total segment
// 	double segment_length = (segLenMin + extra / numSegs) / 2.0;

// 	//return value for the function
// 	shared_ptr<vector<shared_ptr<segment>>> retVec(new vector<shared_ptr<segment>>(numSegs));
// 	//Must initialize the shared_ptrs in the return vector
// 	for (int i = 0; i < numSegs; i++)
// 	{
// 		segment initSeg;
// 		(*retVec)[i] = make_shared<segment>(initSeg);
// 	}

// 	//create a starting position for the segments and set it to first point
// 	Vector3d firstPos = calcEndPt(0, -cylinderHeight/2.0);

// 	//Length that has been covered since the end of the previous section
// 	double currLen = 0;
// 	int currSeg = 0;
// 	bool finalSeg = false;
// 	//double segLenDeb = 0; //Debug parameter for total segment length checking
// 	//calculate the rest of the points for the remaining segments
// 	for (int i = 0; i < number_of_points && currSeg < numSegs; i++)
// 	{
// 		//need to initialize the tbl vector otherwise nothing can be assigned to it
// 		((*retVec)[currSeg])->tbl = make_shared<vector<tableElem>>(vector<tableElem>(0));
// 		((*retVec)[currSeg])->rateVec = make_shared<vector<double>>(vector<double>(0));
// 		((*retVec)[currSeg])->segNum = currSeg;
		
// 		((*retVec)[currSeg])->p1 = firstPos;
// 		Vector3d currPos = firstPos; //curr point being analyzed
// 		Vector3d nextPos; //next curr point to be analyzed
// 		double currSecLen = 0; //amount of space between currPos and nextPos
// 		//segLenDeb = 0;
// 		while (currLen < segment_length)
// 		{
// 			nextPos = getPoint(i); //get next point
// 			currSecLen = (currPos - nextPos).norm();
// 			currLen += currSecLen;//add length due to point to curLen
// 			//segLenDeb += currSecLen;
// 			currPos = getPoint(i); //set up for next iteration
// 			i++; //move to next point
// 		}
// 		i -= 2; //i incremented at end of while and prev inc was too much, so move back 2
// 		currLen -= currSecLen; //take off last addition
// 		//segLenDeb -= currSecLen;
// 		extra = segment_length - currLen;
// 		if (extra < 0)
// 		{
// 			cout << "Calculation of segments failed due to negative extra parameter!!!" << endl;
// 			exit(EXIT_FAILURE);
// 		}
// 		currPos = calcEndPt(i, extra);
// 		//segLenDeb += (currPos - getPoint(i)).norm();
// 		((*retVec)[currSeg])->mid = currPos;
// 		i++;

// 		currLen = 0; //reset the current seg length
// 		while (currLen < segment_length)
// 		{
// 			nextPos = getPoint(i); //get next point
// 			currSecLen = (currPos - nextPos).norm();
// 			currLen += currSecLen;//add length due to point to curLen
// 			//segLenDeb += currSecLen;
// 			currPos = getPoint(i); //set up for next iteration
// 			i++; //move to next point
// 			if (i == number_of_points) //checking for final point, ==number_of_points because of i++
// 			{
// 				finalSeg = true;
// 				break;
// 			}
// 		}
// 		//proceed as usual if not the final segment
// 		if (!finalSeg){
// 			i -= 2; //i incremented at end of while and prev inc was too much, so move back 2
// 			currLen -= currSecLen; //take off last addition
// 			//segLenDeb -= currSecLen;
// 			extra = segment_length - currLen;
// 			if (extra < 0)
// 			{
// 				cout << "Calculation of segments failed due to negative extra parameter.\n";
// 				exit(EXIT_FAILURE);
// 			}
// 			firstPos = calcEndPt(i, extra);
// 			//segLenDeb += (firstPos - getPoint(i)).norm();
// 			((*retVec)[currSeg])->p2 = firstPos;
// 			currLen = 0;
// 		} 
// 		else //final segment needs to 
// 		{
// 			((*retVec)[currSeg])->p2 = calcFinalEndPt(i);
// 		}
// 		currSeg++;
// 	}

// 	if (!(currSeg == numSegs))
// 	{
// 		cout << "Error: Number of calculated segments not as expected.\n";
// 		exit(EXIT_FAILURE);
// 	}

// 	return retVec;
// }

// Vector3d CNT::getPoint(int idx)
// {
// 	Vector3d retVec;
// 	if (idx < 0 || idx > number_of_points - 1)
// 	{
// 		cout << "Invalid index used to access CNT position data!!!" << endl;
// 		exit(EXIT_FAILURE);
// 	}
	
// 	retVec(0, 0) = positions[0][idx];
// 	retVec(1, 0) = positions[1][idx];
// 	retVec(2, 0) = positions[2][idx];

// 	return retVec;
// }

// Vector3d CNT::calcEndPt(int idx, double extra)
// {
// 	Vector3d retVec;
// 	if (idx < 0 || idx > number_of_points - 1)
// 	{
// 		cout << "Invalid index used to access CNT position data!!!" << endl;
// 		exit(EXIT_FAILURE);
// 	}

// 	Vector3d r1 = getPoint(idx);
// 	Vector3d slope = getPoint(idx+1) - r1;
// 	//calculations checked and are correct
// 	retVec = r1 + slope*(2*extra / (cylinderHeight + tubeSeparation));
// 	return retVec;
// }

// Vector3d CNT::calcFinalEndPt(int idx)
// {
// 	Vector3d retVec;
// 	if (idx < 0 || idx > number_of_points)
// 	{
// 		cout << "Invalid index used to access CNT position data!!!" << endl;
// 		exit(EXIT_FAILURE);
// 	}

// 	Vector3d r1 = getPoint(idx-1);
// 	Vector3d slope = r1 - getPoint(idx - 2);
// 	//calculations checked and are correct
// 	retVec = r1 + slope*(cylinderHeight / (cylinderHeight + tubeSeparation));
// 	return retVec;
// }

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

	segments = make_shared<vector<shared_ptr<segment>>>();

	double chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *c;
	gsl_multifit_linear_workspace *work;

	gsl_vector *first_point;
	gsl_vector *second_point;


	print_segment_points(0, number_of_points);
	exit(EXIT_SUCCESS);

	double current_length = 0.;

	int first_idx, second_idx;

	first_idx = 0;

	for (int i = 0; i<number_of_points-1; i++)
	{

		gsl_vector *my_vector = gsl_vector_alloc(3);

		for (int dim = 0; dim < 3; dim++)
		{
			// gsl_vector_set(my_vector, dim, positions[dim][i]-positions[dim][i+1]);
			double delta = gsl_matrix_get(coordinates, i, dim) - gsl_matrix_get(coordinates, i+1, dim);
			gsl_vector_set(my_vector, dim, delta);
		}


		double tmp_double = gsl_blas_dnrm2(my_vector);
		if (tmp_double > segment_length)
		{
			cout << "segment_length is too small for current CNT spine mesh!!!" << endl;
			exit(EXIT_FAILURE);
		}
		current_length += tmp_double;


		if (current_length > segment_length)
		{
			second_idx = i+1;
			int n = second_idx-first_idx;

			X = gsl_matrix_alloc (n, 3);
			y = gsl_vector_alloc (n);
			c = gsl_vector_alloc (3);
			cov = gsl_matrix_alloc (3, 3);
			work = gsl_multifit_linear_alloc (n, 3);

			for (int j = 0; j<n; j++)
			{

				gsl_matrix_set(X, j, 0, 1.0);
				gsl_matrix_set(X, j, 1, gsl_matrix_get(coordinates, j+first_idx, 0));
				gsl_matrix_set(X, j, 2, gsl_matrix_get(coordinates, j+first_idx, 1));

				gsl_vector_set(y, j, gsl_matrix_get(coordinates, j+first_idx, 2));
			}
			gsl_multifit_linear(X, y, c, cov, &chisq, work);

			print_segment_points(first_idx, n);

			exit(EXIT_SUCCESS);

			gsl_matrix_free(X);
			gsl_vector_free(y);
			gsl_vector_free(c);
			gsl_matrix_free(cov);
			gsl_multifit_linear_free(work);

			first_idx = second_idx;
			current_length = 0.;
		}

		gsl_vector_free (my_vector);
	}

	stringstream log_input;
	log_input << "current_length = " << std::scientific << current_length << "   segment_length = " << std::scientific << segment_length <<"    length = " << std::scientific << length;
	write_log(log_input.str());

	exit(EXIT_SUCCESS);

}

gsl_vector* CNT::get_position(int n)
{
	gsl_vector* position_vector = gsl_vector_alloc(3);

	for (int dim = 0; dim < 3; dim++)
	{
		gsl_vector_set(position_vector, dim, positions[dim][n]);
	}

	return position_vector;

}


void CNT::print_segment_points(int first_idx, int n)
{

	stringstream log_input;

	// log_input << "first_idx = " << first_idx << "    length = " << n << endl;

	for (int j = 0; j<n; j++)
	{
		log_input 	<< std::scientific << gsl_matrix_get(coordinates, j+first_idx, 0) << "    "
					<< std::scientific << gsl_matrix_get(coordinates, j+first_idx, 1) << "    "
					<< std::scientific << gsl_matrix_get(coordinates, j+first_idx, 2)
					<< endl;
	}

	// write_log(log_input.str());

	ofstream my_file;
	my_file.open("bullet_physics_points.dat", ios::trunc);
	my_file << log_input.str();

	if (my_file.fail())
	{
		cout << "error in writing to a file!!!" << endl;
		exit(EXIT_FAILURE);
	}

	my_file.close();
}