/**
CNT.cpp
Stores all relevant information for a carbon nanotube

@author Alex Gabourie
@version 1.00
*/

#include "stdafx.h"
#include "CNT.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <regex>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

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
	//positions already has default value

}
/**
Reads a CNT file and creates a CNT object with all the information stored
in that file.

@param filePath The path of the file containing the CNT info
@return CNT Object
*/
CNT::CNT(const string fileName, const string folderPath)
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
			cout << "Error: Incorrect file names. File names should look like \"CNT_Num_x.csv\""
				" where\n \'x\' is a decimal number.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}
	

	//Open .csv file to read.
	ifstream file(filePath);
	//Checks if file is open and readable
	if (!file.good())
	{
		cout << "Cannot read " + fileName + "\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	//Can now read the file
	string temp = " "; //stores intermediate strings during parsing

	//Chirality, first line//
	getline(file, temp, '\n');
	{
		istringstream ss(temp);
		string chir = " ";
		getline(ss, chir, ',');
		getline(ss, chir, ',');
		rgx.assign("(\\d+)(\\s+)(\\d+)");
		smatch matches; //match_results for string objects
		regex_match(chir, matches, rgx);
		if (!matches.empty())
		{
			try{
				n = stoi(matches[1]);
				m = stoi(matches[3]);
			} 
			catch (runtime_error err)
			{
				cout << err.what();
				cout << "\n";
				system("pause");
				exit(EXIT_FAILURE);
			}
		}
		//cannot extract chirality 
		else
		{
			cout << "Error: Cannot extract chirality.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}

	//CNT Length, second line
	getline(file, temp, '\n');
	{
		istringstream ss(temp);
		string len_string = " ";
		getline(ss, len_string, ',');
		getline(ss, len_string, ',');
		try
		{
			length = stod(len_string);
		}
		catch (runtime_error err)
		{
			cout << err.what();
			cout << "\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}

	//Cylinder height, third line
	getline(file, temp, '\n');
	{
		istringstream ss(temp);
		string cyl_string = " ";
		getline(ss, cyl_string, ',');
		getline(ss, cyl_string, ',');
		try
		{
			cylinderHeight = stod(cyl_string);
		}
		catch (runtime_error err)
		{
			cout << err.what();
			cout << "\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}

	//Intertube spacing, fourth line
	getline(file, temp, '\n');
	{
		istringstream ss(temp);
		string tube_space_string = " ";
		getline(ss, tube_space_string, ',');
		getline(ss, tube_space_string, ',');
		try
		{
			minSpacing = stod(tube_space_string);
		}
		catch (runtime_error err)
		{
			cout << err.what();
			cout << "\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}
	
	//intercylinder spacing, fifth line
	getline(file, temp, '\n');
	{
		istringstream ss(temp);
		string cyl_space_string = " ";
		getline(ss, cyl_space_string, ',');
		getline(ss, cyl_space_string, ',');
		try
		{
			tubeSeparation = stod(cyl_space_string);
		}
		catch (runtime_error err)
		{
			cout << err.what();
			cout << "\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}
	
	//Now that all parameters are extracted, calculate diameter
	setDiameter(n, m);

	//move past coordinate labels line
	getline(file, temp, '\n');

	//Counting line loop
	int posNum = -1; //count of positions, negative b/c while loop exits
						//after going an extra line
	while (!file.eof())
	{
		getline(file, temp,'\n');
		posNum++;
	} 

	//With posNum defined, can initialize the positions array
	for (int i = 0; i < 3; i++)
	{
		positions[i].resize(posNum);
	}

	//Position arrays have been defined. Restart file and go to positions again
	file.clear(); //reset state flags
	file.seekg(0); //move to beginning of file
	//move to correct line to start position parsing
	for (int i = 0; i < 6; i++)
	{
		getline(file, temp, '\n');
	}

	string currPos = " ";
	for (int i = 0; i < posNum && file.good(); i++)
	{
		//get next line
		getline(file, temp, '\n');
		//check if we're done
		if (file.eof())
		{
			break;
		}
		
		//Grab all of the position data
		try
		{
			istringstream ss(temp);
			getline(ss, currPos, ',');
			positions[0][i] = stod(currPos);
			getline(ss, currPos, ',');
			positions[1][i] = stod(currPos);
			getline(ss, currPos, ',');
			positions[2][i] = stod(currPos);
		} 
		catch (runtime_error err)
		{
			cout << err.what();
			cout << "\n";
			system("pause");
			exit(EXIT_FAILURE);
		}

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
 