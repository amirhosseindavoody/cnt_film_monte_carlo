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
	numPt = 0;
	positions = vector<vector<double>>(3);

}
/**
Reads a CNT file and creates a CNT object with all the information stored
in that file.

@param filePath The path of the file containing the CNT info
@return CNT Object
*/
CNT::CNT(const string fileName, const string folderPath, double segLen)
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
	numPt = posNum; //assign count to instance variable
	//With posNum defined, can initialize the positions array
	positions = vector<vector<double>>(3);
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
	//Calculate the segments needed for table generation
	segs = calculateSegments(segLen);

	//Checks to see if some segments were calculated
	if (segs->empty())
	{
		cout << "Error: No segments calculated for tube number: ";
		cout << cntNum;
		cout << "\n";
		system("pause");
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

@param segLen The desired length of the segments
@return Vector of the segments
*/
shared_ptr<vector<segment>> CNT::calculateSegments(double segLenMin)
{
	//parameter check
	if (length < segLenMin)
	{
		cout << "Error: Tube length is smaller than minimum segment length.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	else if (segLenMin <= 0)
	{
		cout << "Error: Minimum segment length must be positive.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	//number of segments to use
	int numSegs = static_cast<int>( length / segLenMin); 
	//extra length past numSegs*segLenMin, important for expanding segLen
	double extra = length - segLenMin*numSegs;
	//The equally lengthed segment lengths, half the length of total segment
	double segLen = (segLenMin + extra / numSegs) / 2.0;

	//return value for the function
	shared_ptr<vector<segment>> retVec(new vector<segment>(numSegs));

	//create a starting position for the segments and set it to first point
	Vector3d firstPos = calcEndPt(0, -cylinderHeight/2.0);

	//Length that has been covered since the end of the previous section
	double currLen = 0;
	int currSeg = 0;
	bool finalSeg = false;
	//double segLenDeb = 0; //Debug parameter for total segment length checking
	//calculate the rest of the points for the remaining segments
	for (int i = 0; i < numPt && currSeg < numSegs; i++)
	{
		//need to initialize the tbl vector otherwise nothing can be assigned to it
		(*retVec)[currSeg].tbl = make_shared<vector<tableElem>>(vector<tableElem>(0));
		(*retVec)[currSeg].segNum = currSeg;
		
		(*retVec)[currSeg].p1 = firstPos;
		Vector3d currPos = firstPos; //curr point being analyzed
		Vector3d nextPos; //next curr point to be analyzed
		double currSecLen = 0; //amount of space between currPos and nextPos
		//segLenDeb = 0;
		while (currLen < segLen)
		{
			nextPos = getPoint(i); //get next point
			currSecLen = (currPos - nextPos).norm();
			currLen += currSecLen;//add length due to point to curLen
			//segLenDeb += currSecLen;
			currPos = getPoint(i); //set up for next iteration
			i++; //move to next point
		}
		i -= 2; //i incremented at end of while and prev inc was too much, so move back 2
		currLen -= currSecLen; //take off last addition
		//segLenDeb -= currSecLen;
		extra = segLen - currLen;
		if (extra < 0)
		{
			cout << "Calculation of segments failed due to negative extra parameter.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
		currPos = calcEndPt(i, extra);
		//segLenDeb += (currPos - getPoint(i)).norm();
		(*retVec)[currSeg].mid = currPos;
		i++;

		currLen = 0; //reset the current seg length
		while (currLen < segLen)
		{
			nextPos = getPoint(i); //get next point
			currSecLen = (currPos - nextPos).norm();
			currLen += currSecLen;//add length due to point to curLen
			//segLenDeb += currSecLen;
			currPos = getPoint(i); //set up for next iteration
			i++; //move to next point
			if (i == numPt) //checking for final point, ==numPt because of i++
			{
				finalSeg = true;
				break;
			}
		}
		//proceed as usual if not the final segment
		if (!finalSeg){
			i -= 2; //i incremented at end of while and prev inc was too much, so move back 2
			currLen -= currSecLen; //take off last addition
			//segLenDeb -= currSecLen;
			extra = segLen - currLen;
			if (extra < 0)
			{
				cout << "Calculation of segments failed due to negative extra parameter.\n";
				system("pause");
				exit(EXIT_FAILURE);
			}
			firstPos = calcEndPt(i, extra);
			//segLenDeb += (firstPos - getPoint(i)).norm();
			(*retVec)[currSeg].p2 = firstPos;
			currLen = 0;
		} 
		else //final segment needs to 
		{
			(*retVec)[currSeg].p2 = calcFinalEndPt(i);
		}
		currSeg++;
	}

	if (!(currSeg == numSegs))
	{
		cout << "Error: Number of calculated segments not as expected.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	return retVec;
}

Vector3d CNT::getPoint(int idx)
{
	Vector3d retVec;
	if (idx < 0 || idx > numPt - 1)
	{
		cout << "Invalid index used to access CNT position data.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	
	retVec(0, 0) = positions[0][idx];
	retVec(1, 0) = positions[1][idx];
	retVec(2, 0) = positions[2][idx];

	return retVec;
}

Vector3d CNT::calcEndPt(int idx, double extra)
{
	Vector3d retVec;
	if (idx < 0 || idx > numPt - 1)
	{
		cout << "Invalid index used to access CNT position data.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	Vector3d r1 = getPoint(idx);
	Vector3d slope = getPoint(idx+1) - r1;
	//calculations checked and are correct
	retVec = r1 + slope*(2*extra / (cylinderHeight + tubeSeparation));
	return retVec;
}

Vector3d CNT::calcFinalEndPt(int idx)
{
	Vector3d retVec;
	if (idx < 0 || idx > numPt)
	{
		cout << "Invalid index used to access CNT position data.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	Vector3d r1 = getPoint(idx-1);
	Vector3d slope = r1 - getPoint(idx - 2);
	//calculations checked and are correct
	retVec = r1 + slope*(cylinderHeight / (cylinderHeight + tubeSeparation));
	return retVec;
}