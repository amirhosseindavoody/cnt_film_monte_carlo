// CNTMC_Project.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <string>
#include <sys/stat.h>
#include "dirent.h"
#include <list>
#include "CNT.h"
#include <vector>
#include <memory>
#include "tableElem.h"
#include <locale>
#include <math.h>
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"


using namespace std;

#define TFAC 2.302585092994046 //Factor to multiply 1/gamma by to get 90% of rand numbers giving tr less than deltaT

//method declarations
string folderPathPrompt(bool incorrect);
string xmlFilePathPrompt(bool incorrect);
string outputFolderPathPrompt(bool incorrect);
string checkPath(string path, bool folder);
double updateSegTable(shared_ptr<vector<CNT>> CNT_List, vector<shared_ptr<segment>>::iterator seg, 
	double maxDist, shared_ptr<vector<vector<int>>> colorMap, shared_ptr<vector<double>> rs, shared_ptr<vector<double>> thetas);
int getIndex(shared_ptr<vector<double>> vec, double val);
int getIndex(shared_ptr<vector<double>> vec, double val, int left, int right);
double getRand(bool excludeZero);
void addSelfScattering(shared_ptr<vector<CNT>> CNT_List, double maxGam);
void assignNextState(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> e, double gamma);
double convertUnits(string unit, double val);
shared_ptr<vector<double>> linspace(double low, double high, int num);
void initRandomNumGen();

int main(int argc, char *argv[])
{
	//Initialize random number generator before anything to ensure that getRand() always works
	initRandomNumGen();

	bool done = false; //Reused boolean variable for looping
	string resultFolderPath = " ";
	string inputXMLPath = " ";
	string outputFolderPath = " ";

	if (argc == 1)
	{
		resultFolderPath = folderPathPrompt(false);
		inputXMLPath = xmlFilePathPrompt(false);
		outputFolderPath = outputFolderPathPrompt(false);
	}
	else if (argc == 2)
	{
		resultFolderPath = argv[1];
		inputXMLPath = xmlFilePathPrompt(false);
		outputFolderPath = outputFolderPathPrompt(false);
	}
	else if (argc == 3)
	{
		resultFolderPath = argv[1];
		inputXMLPath = argv[2];
		outputFolderPath = outputFolderPathPrompt(false);
	}
	else if (argc > 4)
	{
		cout << "Incorrect parameters. Only enter file path of results folder to be processed.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		resultFolderPath = argv[1];
		inputXMLPath = argv[2];
		outputFolderPath = argv[3];
	}

	resultFolderPath = checkPath(resultFolderPath, true); //check result folder path
	inputXMLPath = checkPath(inputXMLPath, false); //check xml file path

	////////////////////////////////// OUTPUT FOLDERS ///////////////////////////////

	string timeStamp = "/Date_";
	string response;
	string outputPath;
	{
		time_t timer;
		struct tm currTime;
		if (time(&timer) != -1)
		{
			errno_t err = localtime_s(&currTime, &timer);
			if (err)
			{
				printf("Invalid argument to localtime.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
		}

		timeStamp = timeStamp + to_string(currTime.tm_mday) + "." + to_string(currTime.tm_mon + 1) +
			"." + to_string(currTime.tm_year % 100) + "_Time_" + to_string(currTime.tm_hour) + "." +
			to_string(currTime.tm_min) + "." + to_string(currTime.tm_sec) + "/";
	}

	outputPath = outputFolderPath + timeStamp;
	wstring wide_string(outputPath.begin(), outputPath.end());
	if (CreateDirectory(wide_string.c_str(), nullptr) == 0)
	{
		auto error = GetLastError();
		if (error == ERROR_ALREADY_EXISTS)
		{
			printf("Output folder already exists. Continue anyways? [y/n]\n");
			cin >> response;
			if (response.compare(string("y")) != 0)
			{
				exit(EXIT_FAILURE);
			}
		}
		else if (error == ERROR_PATH_NOT_FOUND)
		{
			printf("Invalid output folder path.\n");
			system("pause");
			exit(EXIT_FAILURE);
		}
	}

	//Results folder exists and can be accessed

	//grab file list
	DIR *resDir;
	struct dirent *ent = nullptr;
	unique_ptr<list<string>>  fileList(new list<string>(0));

	/*
	The next block creates a list of files in the directory chosen to be looked at.
	*/

	//////////////////////////// BUILD FILE LIST ///////////////////////////////////////////

	//Check if folder can be opened - should work due to above checks
	if ((resDir = opendir(resultFolderPath.c_str())) != nullptr)
	{
		//throw away first two results as they are . and ..
		readdir(resDir); readdir(resDir);
		//iterate over all of the real files
		while ((ent = readdir(resDir)) != nullptr)
		{
			fileList->push_back(ent->d_name);
		}
		closedir(resDir); //deletes pointer
	}
	else
	{
		cout << "Could not open directory. Please try program again.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	delete ent;

	/////////////////////////////////// XML FILE PARSE /////////////////////////////////////
	double rmax = 0; //maximum possible differences between sections of CNTs
	while (!done)
	{
		try
		{
			rapidxml::xml_document<> doc; //create xml object
			rapidxml::file<> xmlFile(inputXMLPath.c_str()); //open file
			doc.parse<0>(xmlFile.data()); //parse contents of file
			rapidxml::xml_node<>* currNode = doc.first_node(); //gets the node "Document" or the root node
			currNode = currNode->first_node(); //Output folder
			//Speed up to device dimensions
			for (int i = 0; i < 6; i++)
			{
				currNode = currNode->next_sibling();
			}
			// DEVICE DIMENSIONS NODE //
			auto xdim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value()));
			auto ydim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->value()));
			auto zdim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->next_sibling()->value()));
			//incorrect units
			if (xdim == INT_MIN / 2.0 || ydim == INT_MIN / 2.0 || zdim == INT_MIN / 2.0)
			{
				printf("Configuration Error: Incorrect units for device dimensions.\nRefer to manual"
					" for valid unit entries.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			//incorrect range
			else if (xdim <= 0 || ydim <= 0 || zdim <= 0)
			{
				printf("Configuration Error: Must enter positive device dimensions.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END DEVICE DIMENSIONS NODE //

			rmax = max(max(xdim, ydim), zdim);

			// Will not let me delete, says it is null pointer even though I can get values from it. 
			//delete currNode; 
			done = true;
		}
		catch (runtime_error err)
		{
			string temp;
			int xmlArrayLength = 260; //Maximum path length for Windows
			cout << err.what();
			cout << "\n";
			cout << "Continue? [y/n]: ";
			cin >> temp;
			if (temp.compare("y") != 0)
			{
				system("pause");
				exit(EXIT_FAILURE);
			}
			char *inputXMLPathArray = new char[xmlArrayLength];
			system("cls");
			cout << "Enter config xml path (Example in program files directory):\n";
			cin.ignore();
			cin.getline(inputXMLPathArray, xmlArrayLength);
			inputXMLPath = inputXMLPathArray;
			delete inputXMLPathArray;
		}
	}

	//////////////////////////// CREATE DIST VECTORS AND ARRAYS ///////////////////////////////

	shared_ptr<vector<double>> rs; //Vector containing range of r's
	rmax = 100 * ceil(rmax / 100.0);
	int numBins = static_cast<int>(rmax / 10.0); //number of bins to place r's into 
	double minBin = rmax / static_cast<double>(numBins); //[Angstroms] The size of the bins
	rs = linspace(minBin, rmax, numBins); //Builds rs vector within valid r range

	//builds angle vector from 1 to 90 degrees. Enough bins to cover all relevant angles
	// do not forget to use radians
	double lowAng = 1 * M_PI / 180.0;
	double highAng = 90 * M_PI / 180.0;
	int numAng = 90; //Number of angles to record. One per degree
	shared_ptr<vector<double>> thetas = linspace(lowAng, highAng, numAng);

	shared_ptr<vector<vector<int>>> heatMap(new vector<vector<int>>(rs->size()));
	//initialize all other vectors in the vector to the correct size
	for (vector<vector<int>>::iterator it = heatMap->begin(); it != heatMap->end(); ++it)
	{
		it->resize(thetas->size());
	}

	/*
	This section creates the list of CNTs with all of the relevant information provided in their
	respective files getting processed.
	*/

	//////////////////////////// BUILD CNT AND SEGS //////////////////////////////////////////

	//Iterate through the files and extract
	shared_ptr<vector<CNT>> CNT_List(new vector<CNT>(0));

	for (list<string>::iterator it = fileList->begin(); it != fileList->end(); ++it)
	{
		CNT_List->push_back(CNT(*(it), resultFolderPath, 100.0));
	}

	//Extra check to ensure that all initilizations were successful
	for (vector<CNT>::iterator it = CNT_List->begin(); it != CNT_List->end(); ++it)
	{
		if (!(*it).isInitialized())
		{
			cout << "CNT initialization failure.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}

	/*
	Building the table: Each CNT has a list of segments which has an empty list of table element objects.
	The next section will be filling the segments vector of table elements by calculating the necessary
	additions to it.
	*/

	//////////////////////////// BUILD TABLE ////////////////////////////////////////////////

	//iterate through all of the CNTs and segments
	double maxDist = 500; //[Angstroms]
	//The total number of segments in the simulation, used in exciton placement
	auto numSegs = 0;
	//The maximum of sums of gammas from each segment. This sets the constant gamma value for the entire simulation
	double gamma = 0;
	//loop over CNTs
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		double newGamma;
		//loop over segments in each CNTs
		for (vector<shared_ptr<segment>>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			//get add to each segment relevant table entries
			newGamma = updateSegTable(CNT_List, segit, maxDist, heatMap, rs, thetas);
			if (newGamma > gamma){ gamma = newGamma; }
			numSegs++;
		}
	}

	/////////////////////////////// OUTPUT HEATMAP //////////////////////////////////////

	string fileName = outputPath + "heatMap.csv";
	ofstream file;
	file.open(fileName);
	file << "R Linspace:," << minBin << "," << rmax << "," << numBins << "\n";
	file << "Theta Linspace:," << lowAng << "," << highAng << "," << numAng << "\n";
	file << "Map (r vs. theta):\n";
	//iterate through all of the rs and then thetas while printing to file
	for (int i = 0; i < rs->size(); i++)
	{
		for (int j = 0; j < thetas->size(); j++)
		{
			file << (*heatMap)[i][j] << ",";
		}
		file << "\n";
	}
	file.close();



	/////////////////////////////// ADD SELF SCATTERING //////////////////////////////////

	addSelfScattering(CNT_List, gamma);

	/*
	Now that the tables have been built, the next step is to populate the mesh with excitons. The way this will happen
	is, after a specific number of excitons are chosen to be created, each exciton will be created and assigned an index
	that corresponds to the CNT mesh. First a random CNT will be chosen. This choice will come from a weighted probability
	distribution, the weight being placed on the number of segments each CNT has vs the total number of segments. After the
	CNT number is chosen, then the segment will be chosen from the segments part of the CNT.
	*/

	//////////////////////////// PLACE EXCITON RANDOMLY //////////////////////////////////////////

	//Build the probability vector for CNTs
	shared_ptr<vector<double>> cntProb(new vector<double>(CNT_List->size()));

	//Must set first index so that rest of the table can be built with simple loop
	vector<CNT>::iterator it = CNT_List->begin();
	(*cntProb)[0] = it->segs->size() / static_cast<double>(numSegs);
	++it;

	auto cntIdx = 1; //keeps index of cntProb vector as iterate through CNT_List
	for ( it; it != CNT_List->end(); ++it)
	{
		(*cntProb)[cntIdx] = it->segs->size() / static_cast<double>(numSegs)+(*cntProb)[cntIdx - 1];
		cntIdx++;
	}

	int numExcitons = 10; //number of excitons to add to simulation

	if (numExcitons > 2*numSegs)
	{
		cout << "Error: More excitons than possible number of places for them in simulation.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	
	//vector that stores all excitons
	shared_ptr<vector<shared_ptr<exciton>>> excitons(new vector<shared_ptr<exciton>>(numExcitons)); 

	//assigns excitons to initial locations
	for (UINT32 exNum = 0; exNum < excitons->size(); exNum++)
	{
		done = false; // Flag for successful exciton assignment
		(*excitons)[exNum] = make_shared<exciton>(exciton());
		//go until suitable position for exciton is found
		while (!done)
		{
			//randomly sets energy level to 1 or 2
			(*excitons)[exNum]->setEnergy(static_cast<int>(round(getRand(false)) + 1));
			(*excitons)[exNum]->setCNTidx(getIndex(cntProb, getRand(false)));
			//sets seg index to a number between 0 and the number of segs for the cnt - 1.
			//  The complicated part gets the number of segments from the current cnt
			(*excitons)[exNum]->setSegidx(rand() % ((*CNT_List)[(*excitons)[exNum]->getCNTidx()].segs->size()));
			shared_ptr<vector<shared_ptr<segment>>> temp((*CNT_List)[(*excitons)[exNum]->getCNTidx()].segs); //get segment list
			done = ((*temp)[(*excitons)[exNum]->getSegidx()])->setExciton((*excitons)[exNum]);
		}
	}

	//////////////////////////////////// TIME STEPS ///////////////////////////////////
	double deltaT = (1 / gamma)*TFAC; //time steps at which statistics are calculated
	double Tmax = deltaT * 1000; //maximum simulation time
	double T = 0; //Current simulation time, also time at which next stats will be calculated

	/*
	This section will consist of iterating until the maximum time has been reached. Each iteration
	for T will contain a single/multiple step for each exciton. 
	*/
	while (T <= Tmax)
	{
		T += deltaT; //set new time checkpoint
		for (UINT32 exNum = 0; exNum < excitons->size(); exNum++)
		{
			double tr_tot = 0; //the sum of all tr's in the current deltaT time step
			while (tr_tot <= deltaT)
			{
				tr_tot+= -(1 / gamma)*log(getRand(true)); // add the tr calculation to current time for individual particle
				if (tr_tot > deltaT)
				{
					//Do data recording
				}
				//choose new state
				assignNextState(CNT_List, (*excitons)[exNum],gamma);
			}
		}
	}
	return 0;
}

/**
Assigns the specified exciton to the next state in the simulation.

@param CNT_List The list of carbon nanotubes
@param e The exciton to be updated
*/
void assignNextState(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> e, double gamma)
{
	//Segment the current exciton is located on
	shared_ptr<segment> seg = (*((*CNT_List)[e->getCNTidx()].segs))[e->getSegidx()]; 
	//Get the table index for the
	int tblIdx = getIndex(seg->rateVec, getRand(false)*gamma);
	//stores information about the excitons destination
	tableElem tbl = (*seg->tbl)[tblIdx];
	/*
	It was decided that there are no limits on the number of excitons that
	can be on a segment. 7/20/15
	*/
	e->setCNTidx(tbl.getTubeidx());
	e->setSegidx(tbl.getSegidx());
}

/**
Adds to each segments' rate vector the self scattering component of the simulation

@param CNT_List The list of carbon nanotubes
@param maxGam For all segments, the maximum of the sum of rates
*/
void addSelfScattering(shared_ptr<vector<CNT>> CNT_List, double maxGam)
{
	//CNT index
	int i = 0;
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		//segment index
		int j = 0; 
		//loop over segments in each CNTs
		for (vector<shared_ptr<segment>>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			double currGam = (*segit)->rateVec->back();
			if (maxGam > currGam)
			{
				(*segit)->tbl->push_back(tableElem(1.0, 0.0, maxGam - currGam, i, j));
				(*segit)->rateVec->push_back(maxGam);
			}
			j++;
		}
		i++;
	}
}


/**
Gets a random number between 0 and 1

@return The random number between 0 and 1
*/
double getRand(bool excludeZero)
{
	if (excludeZero)
	{
		int r;
		while ((r = rand()) == 0){}
		return static_cast<double>(r) / static_cast<double>(RAND_MAX);
	}
	return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

/**
Finds the index of the vector that has the number closest to but greater than val.
Does this recursively.

@param vec The vector to search
@param prob The number to compare the indicies to
@return the index that requires the conditions in the method description
*/
int getIndex(shared_ptr<vector<double>> vec, double val)
{	
	//parameter check
	if (vec == nullptr)
	{
		cout << "Error: Empty vector passed to getIndex()";
		system("pause");
		exit(EXIT_FAILURE);
	}
	//small simple cases that do not work with recursive helper method
	if (vec->size() == 1) {return 0;}
	if (vec->size() == 2)
	{
		if ((*vec)[0] >= val){ return 0; }
		return 1;
	}

	//use recursive helper method
	return getIndex(vec, val, 0, vec->size() - 1);
}

/**
The helper method of the 2 parameter get index method. Simple binary search

@param vec The vector to search
@param prob The number to compare the indecies to
@param right The highest index that is part of the vector
@param lef The lowest index that is part of the vector
@return the index that requires the conditions in the method description. Returns -1 if it fails
*/
int getIndex(shared_ptr<vector<double>> vec, double val, int left, int right)
{
	if (right < left) { return left; }

	if (right == left)
	{
		if ((*vec)[right] < val) { return right + 1; }
		return right; 
	} //base case
	
	int center = static_cast<int>(floor(static_cast<double>(right + left) / 2.0));
	if ((*vec)[center] == val){ return center; } //base case
	
	//recursive cases
	if ((*vec)[center] > val){ return getIndex(vec, val, left, center - 1); }
	if ((*vec)[center] < val){ return getIndex(vec, val, center + 1, right); }
	return -1;

}

/**
Takes a segment and determines which elements should be added to its tables based on distance away
from the segment.

@param CNT_List The list of carbon nanotubes to iterate over
@param seg The segment that we will be adding table elements to
@param maxDist If segments are within maxDist of seg, then they will be added to the table
@param colorMap A count of all rs at particular thetas to get mesh statistics
@return The sum of all the rates calculated for the segment. For transition purposes
*/
double updateSegTable(shared_ptr<vector<CNT>> CNT_List, vector<shared_ptr<segment>>::iterator seg,
	double maxDist, shared_ptr<vector<vector<int>>> heatMap, shared_ptr<vector<double>> rs, shared_ptr<vector<double>> thetas)
{
	double rate = 0;
	//iterate over CNTs
	int i = 0; //CNT index counter
	//originally structured without i and j
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		int j = 0; //segment index counter
		//iterate over all segments considered for seg
		for (vector<shared_ptr<segment>>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			double r = tableElem::calcDist((*seg)->mid, (*segit)->mid);
			auto theta = tableElem::calcThet(seg, segit);

			//Color Map Additions
			if (r != 0)
			{
				(*heatMap)[getIndex(rs, r)][getIndex(thetas, theta)]++;
				//Check if within range
				if (r <= maxDist)
				{
					auto g = 6.4000e+19; //First draft estimate
					(*seg)->tbl->push_back(tableElem(r, theta, g, i, j)); //tbl initialized in CNT::calculateSegments
					(*seg)->rateVec->push_back(rate += ((*seg)->tbl->back()).getRate());//tbl initialized in CNT::calculateSegments
				}
			}	
			j++;
		}
		i++;
	}
	return rate;
}


/**
Gets the path of the folder containing CNT mesh results

@param incorrect Runs special prompt in case that cmd args were incorrect
@return The string containing the results folder path.
*/
string folderPathPrompt(bool incorrect)
{
	if (incorrect)
	{
		string temp = " ";
		cout << "Re-enter results folder path? [y/n]: ";
		cin >> temp;
		//check response, end if decline retry
		if (temp.compare("y") != 0)
		{
			system("pause");
			exit(EXIT_FAILURE);
		}

	}
	
	int inputPathLengthMax = 260; //Maximum path length for Windows
	string returnString;
	char *inputResultsFolderPathArray = new char[inputPathLengthMax];
	cout << "Enter path of results folder:\n";
	if (incorrect)
		cin.ignore(); //if reentering, must ignore the next input
	cin.getline(inputResultsFolderPathArray, inputPathLengthMax);
	returnString = inputResultsFolderPathArray;
	delete[] inputResultsFolderPathArray;
	return returnString;
}

/**
Gets the path of the file used to configure the CNT mesh generation

@param incorrect Runs special prompt in case that cmd args were incorrect
@return The string containing the XML file path.
*/
string xmlFilePathPrompt(bool incorrect)
{
	if (incorrect)
	{
		string temp = " ";
		cout << "Re-enter XML config file path? [y/n]: ";
		cin >> temp;
		//check response, end if decline retry
		if (temp.compare("y") != 0)
		{
			system("pause");
			exit(EXIT_FAILURE);
		}
	}

	int inputPathLengthMax = 260; //Maximum path length for Windows
	string returnString;
	char *inputXMLFilePathArray = new char[inputPathLengthMax];
	cout << "Enter path of XML config file:\n";
	if (incorrect)
		cin.ignore(); //if reentering, must ignore the next input
	cin.getline(inputXMLFilePathArray, inputPathLengthMax);
	returnString = inputXMLFilePathArray;
	delete[] inputXMLFilePathArray;
	return returnString;
}

/**
Gets the path of output folder

@param incorrect Runs special prompt in case that cmd args were incorrect
@return The string containing the XML file path.
*/
string outputFolderPathPrompt(bool incorrect)
{
	if (incorrect)
	{
		string temp = " ";
		cout << "Re-enter output folder path? [y/n]: ";
		cin >> temp;
		//check response, end if decline retry
		if (temp.compare("y") != 0)
		{
			system("pause");
			exit(EXIT_FAILURE);
		}
	}

	int inputPathLengthMax = 260; //Maximum path length for Windows
	string returnString;
	char *inputXMLFilePathArray = new char[inputPathLengthMax];
	cout << "Enter path of output folder:\n";
	if (incorrect)
		cin.ignore(); //if reentering, must ignore the next input
	cin.getline(inputXMLFilePathArray, inputPathLengthMax);
	returnString = inputXMLFilePathArray;
	delete[] inputXMLFilePathArray;
	return returnString;
}

/**
Converts numbers with some units to angstroms

@param unit The current unit
@param val The current value
@return the value in angstroms
*/
double convertUnits(string unit, double val)
{
	if (unit.compare("mm") == 0 || unit.compare("millimeter") == 0)
	{
		return val * 10000000;
	}
	else if (unit.compare("um") == 0 || unit.compare("micrometer") == 0)
	{
		return val * 10000;
	}
	else if (unit.compare("nm") == 0 || unit.compare("nanometer") == 0)
	{
		return val * 10;
	}
	else if (unit.compare("pm") == 0 || unit.compare("picometer") == 0)
	{
		return val * .01;
	}
	else if (unit.compare("A") == 0 || unit.compare("angstrom") == 0)
	{
		return val;
	}
	else
	{
		return INT_MIN;
	}
}

/**
Checks to see if the file/folder exists and is accessable

@param path The path of the file/folder
@param folder If true, then the path is for a folder, file otherwise
@return The string of the path that works
*/
string checkPath(string path, bool folder)
{
	//Checking validity of folder path
	while (true)
	{
		//Container for file/folder information
		struct stat info;

		if (stat(path.c_str(), &info) != 0)
		{
			printf("Cannot access %s\n", path.c_str());
			//get new path for next itr or quit
			if (folder){ path = folderPathPrompt(true); }
			else { path = xmlFilePathPrompt(true);}
		}
		//path is accessable
		else
		{
			//Check to see if correctly dir or file
			if ((folder && (info.st_mode & S_IFDIR)) || 
				(!folder && (info.st_mode & S_IFREG))) { return path; }
			
			//Incorrect type for the path. Reprompt or quit
			if (folder)
			{
				printf("%s is not a directory\n", path.c_str());
				path = folderPathPrompt(true);
			} 
			else
			{
				printf("%s is not a file\n", path.c_str());
				path = xmlFilePathPrompt(true);
			}
		}
	}
}

/**
Creates a vector with matlab style linspace numbering

@param low The lowest number to be in vector
@param high The highest number to be in vector
@param num The number of points to be in vector
@return A pointer to the resulting vector
*/
shared_ptr<vector<double>> linspace(double low, double high, int num)
{
	shared_ptr<vector<double>> retVec(new vector<double>(num));
	double step = (high - low) / static_cast<double>(num - 1);
	for (int i = 0; i < num; i++)
	{
		(*retVec)[i] = low;
		low += step;
	}
	return retVec;
}

/**
Initialized the random number generator by providing a seed value from 
the time.
*/
void initRandomNumGen()
{
	//Initialize random number generation
	time_t seconds;
	time(&seconds); //assign time from clock
	//Seed the random number generator
	srand(static_cast<int>(seconds));
}