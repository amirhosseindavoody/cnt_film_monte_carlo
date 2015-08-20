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
#include <regex>
#include <Windows.h>
#include <omp.h>
#include <conio.h>
#include "typeTransition.h"
#include "chirality.h"
#include "paramStructs.h"


using namespace std;

#define MAX_T_STEPS 10000000 //maximum number of time steps allowed for the simulation

//function pointer for a function that takes current data and translates it into a value for the rate table
//Typedef so that when passing as a parameter, the entire declaration doesn't need to be typed in.
typedef void(*dat2tab) (tableUpdater&);

//method declarations
string folderPathPrompt(bool incorrect);
string xmlFilePathPrompt(bool incorrect);
string outputFolderPathPrompt(bool incorrect);
string checkPath(string path, bool folder);
double getRand(bool excludeZero);
double convertUnits(string unit, double val);
shared_ptr<vector<double>> linspace(double low, double high, int num);
void initRandomNumGen();
void injectExciton(shared_ptr<exciton> exciton, shared_ptr<vector<shared_ptr<segment>>> inContact);
void writeStateToFile(shared_ptr<ofstream> file, shared_ptr<vector<int>> currCount, double T);
double diffclock(clock_t end, clock_t start);
string getRunStatus(double T, double Tmax, double runtime, boolean runtimeKnown);
void ClearScreen();
string getRunTime(double runtime);
string fixPath(string &path);
string GetLastErrorAsString();
int getIndex(shared_ptr<vector<double>> vec, double val);
int getIndex(vector<Chirality> &vec, Chirality &val);
void addChiralitiesToList(vector<Chirality> &vec, char* filename);
chirPair getChiralityFromFilename(char* filename);


double updateSegTable(double maxDist, dat2tab addDataToTable, tableUpdater &t, heatMapInfo &h);

void addSelfScattering(shared_ptr<vector<CNT>> CNT_List, double maxGam);
void assignNextState(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> e, double gamma, shared_ptr<vector<double>> regionBdr);
bool hasMovedToOutContact(shared_ptr<exciton> exciton, shared_ptr<vector<double>> regionBdr, shared_ptr<vector<CNT>> CNT_List);
void markCurrentExcitonPosition(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> exciton, shared_ptr<vector<int>> currCount,
	shared_ptr<vector<double>> regionBdr);
void updateExcitonList(int numExcitonsAtCont, shared_ptr<vector<shared_ptr<exciton>>> excitons, shared_ptr<vector<int>> currCount,
	shared_ptr<vector<shared_ptr<segment>>> inContact);
void writeExcitonDistSupportingInfo(bool tableFromFile, string outputPath, int numExcitons, double Tmax, double deltaT, double segLenMin, int numRegions,
	double xdim, double minBin, double rmax, int numBins, double lowAng, double highAng, int numAng, UINT64 numTSteps, double regLenMin, string runtime);
void addDataToTableCalc(tableUpdater &t);
void addDataToTableRead(tableUpdater &t);
void sortChiralities(vector<Chirality> &list);


//Global variables
double ymax = 0; //stores maximum height the cylinders of the CNTs are found at. All will be greater than 0.

//Runs the file input, monte carlo, and file output sections of code
int main(int argc, char *argv[])
{

	//timing functions
	clock_t start = clock();

	//Initialize random number generator before anything to ensure that getRand() always works
	initRandomNumGen();
	
	tableUpdater tableParams;//Paramater struct
	unsigned int NUM_THREADS = omp_get_max_threads();
	string status;
	//Varible initialization
	double segLenMin = 100.0; //[Angstroms]
	double regLenMin = segLenMin;
	double runtime;
	int numSteps = 10000; //number of deltaT's the simulation runs if the user specifies
	double maxDist = 300; //[Angstroms] The maximum segment separation to be considered for exciton transfer
	//The number of excitons to be in injection contact at start
	int numExcitonsAtCont = 100000;
	bool autoComplete = false;
	double threshold = .01; //The value that the difference/maxDiff must be below to finish
	int numToFinish = 5;// Number of differences/maxDiff  that must be below thresh in a row to finish
	int numToCheck = 1000; //number of time steps to check finish completion
	double tfac = log(.3);
	dat2tab addDataToTable; //function pointer for updating table
	bool tableFromFile = false; //tells simulation to either read table from file or create new table
	vector<Chirality> meshChirList = vector<Chirality>(0);

	bool done = false; //Reused boolean variable for looping
	string resultFolderPath = " ";

	if (argc == 1)
	{
		resultFolderPath = folderPathPrompt(false);
	}
	else if (argc > 2)
	{
		cout << "Incorrect parameters. Only enter file path of results folder to be processed.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	else
	{
		resultFolderPath = argv[1];
	}

	resultFolderPath = fixPath(resultFolderPath);
	resultFolderPath = checkPath(resultFolderPath, true); //check result folder path


	////////////////////////////////// OUTPUT FOLDERS ///////////////////////////////

	//Outputfolder should be the same as the result folder for the mesh
	// All checks are already completed
	string outputPath = resultFolderPath + "/"; 

	//Results folder exists and can be accessed
	//grab file list
	DIR *resDir;
	struct dirent *ent = nullptr;
	unique_ptr<list<string>>  fileList(new list<string>(0));

	/*
	The next block creates a list of files in the directory chosen to be looked at.
	*/

	//////////////////////////// BUILD FILE LIST ///////////////////////////////////////////
	regex rgx("CNT_Num_\\d+\\.csv"); //files we want to look through
	//Check if folder can be opened - should work due to above checks
	if ((resDir = opendir(resultFolderPath.c_str())) != nullptr)
	{
		//throw away first two results as they are . and ..
		readdir(resDir); readdir(resDir);
		//iterate over all of the real files
		while ((ent = readdir(resDir)) != nullptr)
		{
			smatch matches; //match_results for string objects
			regex_search(string(ent->d_name), matches, rgx);
			if (!matches.empty())
			{
				fileList->push_back(ent->d_name);
			}
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
	
	//First need to get correct xml file name
	auto xmlNamePos = resultFolderPath.rfind('/', resultFolderPath.size() - 1) + 1;
	string inputXMLPath = outputPath + resultFolderPath.substr(xmlNamePos, outputPath.size() - xmlNamePos) + ".xml";
	
	double rmax = 0; //maximum possible differences between sections of CNTs
	double xdim = 0; //Dimension in which exciton populations will be monitored
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
			xdim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value()));
			auto ydim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->value()));
			auto zdim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->next_sibling()->value()));
			//incorrect units
			if (xdim == INT_MIN || ydim == INT_MIN || zdim == INT_MIN)
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

			//x and z dim are the ground dimensions and y is height. Need to make sure bottom corner
			// to top other corner is included in r range. ymax is not defined so this is intermediate
			// value. It is changed after CNT initialization.
			rmax = sqrt(pow(xdim,2)+pow(zdim,2));

			// READ CHIRALITIES //
			rapidxml::xml_node<>* chirNode = currNode->next_sibling(); //get to chirality
			chirNode = chirNode->first_node();
			while (chirNode)
			{
				meshChirList.push_back(Chirality(atoi(chirNode->first_node()->value()), 
					atoi(chirNode->first_node()->next_sibling()->value())));
				chirNode = chirNode->next_sibling();
			}
			sortChiralities(meshChirList);
			// END READ CHIRALITIES //

			// USE PREBUILT TABLE NODE //
			currNode = currNode->next_sibling()->next_sibling(); //to prebuilt node
			//check to see state
			if (!string(currNode->value()).compare("true"))
			{
				tableFromFile = true;
			}
			// END USE PREBUILT TABLE NODE //

			// NUMBER OF EXCITONS NODE //
			currNode = currNode->next_sibling(); //to number of excitons
			int exNumTemp = atoi(currNode->value());
			//accept only positive number of excitons otherwise have default value
			if (exNumTemp > 0)
			{
				numExcitonsAtCont = exNumTemp;
			}
			else
			{
				printf("Configuration Error: Must have positive number of excitons.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END NUMBER OF EXCITONS NODE //

			// REGION LENGTH NODE //
			currNode = currNode->next_sibling();
			double regLenTemp = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value()));
			if (regLenTemp >= 0 && regLenTemp < xdim)
			{
				regLenMin = regLenTemp;
			}
			else
			{
				printf("Configuration Error: Region length must be positive number.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END REGION LENGTH NODE //

			// SEGMENT LENGTH NODE //
			currNode = currNode->next_sibling();
			double segLenTemp = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value()));
			if (segLenTemp >= 0)
			{
				segLenMin = segLenTemp;
			}
			else
			{
				printf("Configuration Error: Segment length must be positive number.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END SEGMENT LENGTH NODE //

			//I want it to be the case that if either segLenMin or regLenMin = 0, then they assume
			// the value of the other
			if (segLenMin == 0 && regLenMin != 0)
			{
				segLenMin = regLenMin;
			}
			else if (segLenMin != 0 && regLenMin == 0)
			{
				regLenMin = segLenMin;
			} 
			else if (segLenMin == 0 && regLenMin == 0)
			{
				regLenMin = segLenMin = 100.0;
			}

			// SEGMENT SEPARATION //
			currNode = currNode->next_sibling();
			double segSepTemp = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value()));
			if (segSepTemp > 0)
			{
				maxDist = segSepTemp;
			}
			else
			{
				printf("Configuration Error: Segment separation must be a positive number.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END SEGMENT SEPARATION //

			// NUMBER OF TIME STEPS //
			currNode = currNode->next_sibling();
			int numTStepsTemp = atoi(currNode->value());
			if (numTStepsTemp > 0)
			{
				numSteps = numTStepsTemp;
			}
			else if (numTStepsTemp == 0)
			{
				autoComplete = true;
				numSteps = MAX_T_STEPS;
			}
			else
			{
				printf("Configuration Error: Must have positive number time steps. 0 steps if auto complete is desired.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END NUMBER OF TIME STEPS //


			// PERCENT FREE FLIGHT TIMES ABOVE DELTA T //
			currNode = currNode->next_sibling();
			double perc_tfacTemp = atoi(currNode->value()) / 100.0;
			if (perc_tfacTemp <= 100 && perc_tfacTemp > 0)
			{
				tfac = abs(log(perc_tfacTemp));
			}
			else
			{
				printf("Configuration Error: Percent of free flight times above delta T must be greater than 0 and less than or equal to 100.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END PERCENT FREE FLIGHT TIMES ABOVE DELTA T  //

			// AUTO COMPLETE //
			currNode = currNode->next_sibling()->first_node(); //enabled?
			//if auto-complete enabled
			if (!string(currNode->value()).compare("true"))
			{
				autoComplete |= true;
				numSteps = MAX_T_STEPS;

				//to threshold node
				currNode = currNode->next_sibling();
				double thresholdTemp = atof(currNode->value());
				//check for range issues
				if (threshold > 0.0 && threshold < 1.0)
				{
					threshold = thresholdTemp;
				}
				else
				{
					printf("Configuration Error: Threshold must be value greater than 0 and less than 1.\n");
					system("pause");
					exit(EXIT_FAILURE);
				}

				//to numBelowThresholdNode
				currNode = currNode->next_sibling();
				int numToFinishTemp = atoi(currNode->value());
				if (numToFinishTemp > 0) //input validation
				{
					numToFinish = numToFinishTemp;
				}
				else
				{
					printf("Configuration Error: Number below threshold must be greater than 0.\n");
					system("pause");
					exit(EXIT_FAILURE);
				}

				//to numToAverage
				currNode = currNode->next_sibling();
				int numToCheckTemp = atoi(currNode->value());
				if (numToCheckTemp > 0)
				{
					numToCheck = numToCheckTemp;
				}
				else
				{
					printf("Configuration Error: Number to average must be greater than 0.\n");
					system("pause");
					exit(EXIT_FAILURE);
				}
			}

			// END AUTO COMPLETE // 


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
	//@@@@@@@@ TIME UPDATE @@@@@@@@//
	clock_t end = clock();
	runtime = diffclock(end, start);
	ClearScreen();
	status = getRunStatus(0, 0, runtime, !autoComplete);
	cout << status << endl;
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@//

	/*
	This section reads information from the table files if the option is chosen. While reading
	the binary files, a data structure is initializes such that the tables can be read for 
	transition table creation later.
	*/


	if (tableFromFile)
	{

		////////////////// TRANSITION TABLE PARAMETERS ////////////////////////////////////

		uint32_t numChiralities; //Number of different chiralities included in the simulation
		double r_low;
		double r_high;
		double t_low;
		double t_high;
		uint32_t r_size;  //number of r's the rates have been calculated for
		uint32_t theta_size; //number of theta's the rates have been calculated for
		vector<Chirality> ahChirList = vector<Chirality>(0); //chiralities coming from amirhossein's tables

		//READ TABLE DETAILS FILE

		//path of transition rate tables folder
		string tableFolderPath = outputPath + "transfer_rate_tables/";
		string tableDetailFilePath;

		unique_ptr<list<string>>  tableFileList(new list<string>(0));
		regex chirrgx("(\\d+),(\\d+)_(\\d+),(\\d+)\\.bin"); //files we want to look through
		regex detailrgx("details\\.csv");
		regex numrgx(""); //matching numbers in file name
		//Check if folder can be opened - should work due to above checks
		if ((resDir = opendir(tableFolderPath.c_str())) != nullptr)
		{
			//throw away first two results as they are . and ..
			readdir(resDir); readdir(resDir);
			//iterate over all of the real files
			while ((ent = readdir(resDir)) != nullptr)
			{
				smatch matches; //match_results for string objects
				regex_search(string(ent->d_name), matches, chirrgx);
				//check if we matched the file name needed for data files
				if (!matches.empty())
				{
					tableFileList->push_back(ent->d_name);

					// ADD CHIRALITIES FROM FILENAME TO CHIRALITY LIST //

					addChiralitiesToList(ahChirList, ent->d_name);

				}
				else
				{
					smatch newMatch;
					regex_search(string(ent->d_name), newMatch, detailrgx);
					if (!newMatch.empty())
					{
						//define the details path if file exists
						tableDetailFilePath = tableFolderPath + "details.csv";
					}
				}
			}
			closedir(resDir); //deletes pointer
			sortChiralities(ahChirList);
		}
		else
		{
			ClearScreen();
			cout << "Could not open transition_rate_tables directory. Please try program again.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
		delete ent;
		
		//check to make sure details are available.
		if (tableDetailFilePath.empty())
		{
			ClearScreen();
			printf("Table Conversion Error: No details.csv file provided in transition_rate_tables folder.\n");
			system("pause");
			exit(EXIT_FAILURE);
		}

		//read details file.
		ifstream detFile(tableDetailFilePath);
		//file can be read
		if (!detFile.good())
		{
			ClearScreen();
			cout << "Cannot read " + tableDetailFilePath + "\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
		string temp = " "; //stores intermediate strings during parsing
		//skip first line as its only for user readability
		getline(detFile, temp, '\n');
		getline(detFile, temp, ',');
		numChiralities = atoi(temp.c_str());
		getline(detFile, temp, ',');
		r_low = atoi(temp.c_str());
		getline(detFile, temp, ',');
		r_high = atoi(temp.c_str());
		getline(detFile, temp, ',');
		r_size = atoi(temp.c_str());
		getline(detFile, temp, ',');
		t_low = atoi(temp.c_str());
		getline(detFile, temp, ',');
		t_high = atoi(temp.c_str());
		getline(detFile, temp, ',');
		theta_size = atoi(temp.c_str());

		//Check that all three ways of specifying number of chiralities are matching up.
		if (tableFileList->size() != numChiralities*numChiralities)
		{
			ClearScreen();
			printf("Error: Not enough chirality files provided. "
				"For n chiralities, must have n^2 transition files to capture all chirality combinations.\n");
			system("pause");
			exit(EXIT_FAILURE);
		}
		else if (numChiralities != meshChirList.size())
		{
			ClearScreen();
			printf("Error: Number of chiralities specified does not match the number of chiralities in mesh simulation.\n");
			system("pause");
			exit(EXIT_FAILURE);
		}
		else if (ahChirList.size() > meshChirList.size())
		{
			ClearScreen();
			printf("Error: More chiralities in binary files than used in mesh generation.\n");
			system("pause");
			exit(EXIT_FAILURE);
		}
		else if (ahChirList.size() < meshChirList.size())
		{
			ClearScreen();
			printf("Error: Fewer chiralities in binary files than used in mesh generation.\n");
			system("pause");
			exit(EXIT_FAILURE);
		}

		//Compare the chirality vectors to ensure that they are the same
		{
			bool compareChirVec = true; //sum of abs() of return values for chiralities
			for (int i = 0; i < meshChirList.size(); i++)
			{
				compareChirVec = compareChirVec && (ahChirList[i] == meshChirList[i]);
			}
			if (!compareChirVec)
			{
				ClearScreen();
				printf("Error: Mesh chiralities and table chiralities do not match.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
		}

		//Below is the structure of the c2c object to show from what chir to what chir the transition occurs at
		/*
		C1,C2,C3   from
		|  |  |     |
		C1,C1,C1    to
		C2,C2,C2
		C3,C3,C3
		*/
		//main list to store input table information
		tableParams.c2c = make_shared<vector<vector<typeTransition>>>(vector<vector<typeTransition>>(numChiralities));
		//initialize the rest c2c object
		auto chirItr = ahChirList.begin();
		for (auto it = tableParams.c2c->begin(); it != tableParams.c2c->end(); ++it, ++chirItr)
		{
			for (auto chirItr2 = ahChirList.begin(); chirItr2 != ahChirList.end(); ++chirItr2)
			{
				it->push_back(typeTransition(*chirItr, *chirItr2, r_size, theta_size));
			}
		}
		tableParams.r_vec = linspace(r_low,r_high,r_size);
		tableParams.t_vec = linspace(t_low, t_high, theta_size);

		// Parse files
		int src_idx;
		int dest_idx;
		for (auto itr = tableFileList->begin(); itr != tableFileList->end(); ++itr)
		{
			//prep text for parse
			char filename[13];
			strcpy_s(filename, _countof(filename), itr->c_str()); 
			//get pair information
			chirPair pair = getChiralityFromFilename(filename);
			//get index of table to edit
			src_idx = getIndex(ahChirList, pair.c1);
			dest_idx = getIndex(ahChirList, pair.c2);

			//get the energyTransition vector
			auto Elist = (*tableParams.c2c)[src_idx][dest_idx].getETransList();

			double readRate;
			ifstream infile;
			infile.open(tableFolderPath + *itr, ios::binary | ios::in);


			for (int r_idx = 0; r_idx < r_size; r_idx++)
			{
				for (int t_idx = 0; t_idx < theta_size; t_idx++)
				{
					for (auto e_itr = Elist->begin(); e_itr != Elist->end(); ++e_itr)
					{
						infile.read(reinterpret_cast<char *>(&readRate), sizeof(double));
						e_itr->setTableValue(r_idx, t_idx, readRate);
					}
				}
			}
			infile.close();
		}


		addDataToTable = addDataToTableRead;

		/////////////////// END OF TRANSITION TABLE PARAMETERS ///////////////////////////
	}
	else
	{
		addDataToTable = addDataToTableCalc;
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
		CNT_List->push_back(CNT(*(it), resultFolderPath, segLenMin));
	}

	//Extra check to ensure that all initilizations were successful
	for (vector<CNT>::iterator it = CNT_List->begin(); it != CNT_List->end(); ++it)
	{
		if (!(*it).isInitialized())
		{
			ClearScreen();
			cout << "CNT initialization failure.\n";
			system("pause");
			exit(EXIT_FAILURE);
		}
	}
	tableParams.CNT_List = CNT_List; //add pointer to parameter list
	//Set final rmax value
	rmax = sqrt(pow(rmax, 2) + pow(ymax,2));
	

	//@@@@@@@@ TIME UPDATE @@@@@@@@//
	end = clock();
	runtime = diffclock(end, start);
	ClearScreen();
	status = getRunStatus(0, 0, runtime, !autoComplete);
	cout << status << endl;
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@//

	//////////////////////////// CREATE DIST VECTORS AND ARRAYS ///////////////////////////////
	heatMapInfo heatMap;
	//shared_ptr<vector<double>> rs; //Vector containing range of r's
	rmax = 100 * ceil(rmax / 100.0);
	int numBins = static_cast<int>(rmax / 10.0); //number of bins to place r's into 
	double minBin = rmax / static_cast<double>(numBins); //[Angstroms] The size of the bins
	heatMap.rs = linspace(minBin, rmax, numBins); //Builds rs vector within valid r range

	//builds angle vector from 1 to 90 degrees. Enough bins to cover all relevant angles
	// do not forget to use radians
	double lowAng = 1 * M_PI / 180.0;
	double highAng = 90 * M_PI / 180.0;
	int numAng = 90; //Number of angles to record. One per degree
	heatMap.thetas = linspace(lowAng, highAng, numAng);

	heatMap.map = make_shared<vector<vector<int>>>(vector<vector<int>>(heatMap.rs->size()));
	//initialize all other vectors in the vector to the correct size
	for (vector<vector<int>>::iterator it = heatMap.map->begin(); it != heatMap.map->end(); ++it)
	{
		it->resize(heatMap.thetas->size());
	}

	/*
	Building the table: Each CNT has a list of segments which has an empty list of table element objects.
	The next section will be filling the segments vector of table elements by calculating the necessary
	additions to it.
	*/

	//////////////////////////// BUILD TABLE ////////////////////////////////////////////////

	///////////// Variables for placing excitons in the future /////////////
	int numRegions = static_cast<int>(xdim / regLenMin); //number of regions in the simulation
	double extra = xdim - numRegions*regLenMin; //extra amount * numRegions that should be added to regLenMin
	double regLen = regLenMin + extra / numRegions; //The length of each region.

	//The boundary of the rgions in the x direction
	shared_ptr<vector<double>> regionBdr = linspace(regLen - (xdim / 2), xdim / 2, numRegions);
	shared_ptr<vector<int>> segCountPerReg = make_shared<vector<int>>(vector<int>(numRegions)); //To get dist stats
	//List of segments in the first region, which is used as a input contact
	shared_ptr<vector<shared_ptr<segment>>> inContact(new vector<shared_ptr<segment>>(0));
	
	//iterate through all of the CNTs and segments
	//The total number of segments in the simulation, used in exciton placement
	auto numSegs = 0;
	//The maximum of sums of gammas from each segment. This sets the constant gamma value for the entire simulation
	double gamma = 0;

	int buildTblDecPerc = static_cast<int>(CNT_List->size() / 10.0); //ten percent of total number of tubes
	int buildTblCntr = 0;
	//loop over CNTs
	for (auto cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		double newGamma;
		//loop over segments in each CNTs
		for (auto segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			tableParams.seg = segit; //set passing params
			int regIdx = getIndex(regionBdr, (*segit)->mid(0));
			(*segCountPerReg)[regIdx]++; //increment the count based on where segment is
			if (regIdx == 0){ inContact->push_back(*segit); } //First region is injection contact
			//get add to each segment relevant table entries
			newGamma = updateSegTable(maxDist, addDataToTable, tableParams, heatMap);
			if (newGamma > gamma){ gamma = newGamma; }
			numSegs++;
		}
		buildTblCntr++;
		if (buildTblCntr == buildTblDecPerc)
		{
			//@@@@@@@@ TIME UPDATE @@@@@@@@//
			end = clock();
			runtime = diffclock(end, start);
			ClearScreen();
			status = getRunStatus(0, 0, runtime, !autoComplete);
			cout << status << endl;
			//@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
			buildTblCntr = 0;
		}
		
	}


	/////////////////////////////// OUTPUT SEG COUNT PER REGION //////////////////////////////////////

	string segmentCountFileName = outputPath + "segmentCountPerRegion.csv";
	ofstream segmentCountFile;
	segmentCountFile.open(segmentCountFileName);
	segmentCountFile << (*segCountPerReg)[0];
	//iterate through all of the rs and then thetas while printing to file
	for (UINT32 i = 1; i < segCountPerReg->size(); i++)
	{
		segmentCountFile << "," << (*segCountPerReg)[i];
	}
	segmentCountFile.close();

	/////////////////////////////// OUTPUT HEATMAP //////////////////////////////////////

	string heatMapFileName = outputPath + "heatMap.csv";
	ofstream heatMapFile;
	heatMapFile.open(heatMapFileName);
	//iterate through all of the rs and then thetas while printing to file
	for (UINT32 i = 0; i < heatMap.rs->size(); i++)
	{
		for (UINT32 j = 0; j < heatMap.thetas->size(); j++)
		{
			heatMapFile << (*heatMap.map)[i][j] << ",";
		}
		heatMapFile << "\n";
	}
	heatMapFile.close();



	/////////////////////////////// ADD SELF SCATTERING //////////////////////////////////

	addSelfScattering(CNT_List, gamma);

	/*
	Now that the tables have been built, the next step is to populate the mesh with excitons. The way this will happen
	is, after a specific number of excitons are chosen to be created, each exciton will be created and assigned an index
	that corresponds to the CNT mesh. The available locations for excitons to start will be determined by the number of 
	segments that are within the delta x region specified. The exciton assignment will be equally and randomly distributed 
	among the list of segments in the delta x region.
	*/

	//////////////////////////// PLACE EXCITONS ////////////////////////////////////////////

	
	//Vector of excitons. Positions and energies must still be assigned
	shared_ptr<vector<shared_ptr<exciton>>> excitons(new vector<shared_ptr<exciton>>(numExcitonsAtCont));
	for (int exNum = 0; exNum < numExcitonsAtCont; exNum++)
	{
		(*excitons)[exNum] = make_shared<exciton>(exciton()); //initialize exciton at location in exciton list
		injectExciton((*excitons)[exNum], inContact);
	}


	//////////////////////////////////// TIME STEPS ///////////////////////////////////
	double deltaT = (1 / gamma)*tfac; //time steps at which statistics are calculated
	double Tmax = deltaT * numSteps; //maximum simulation time
	double T = 0; //Current simulation time, also time at which next stats will be calculated

	shared_ptr<vector<int>> currCount; //The count of excitons in each region of the CNT mesh

	//File output initializations
	string excitonDistFileName;
	//choose file name based on type of simulation
	if (tableFromFile)
	{
		excitonDistFileName = outputPath + "excitonDist_t.csv";
	}
	else
	{
		excitonDistFileName = outputPath + "excitonDist.csv";
	}
	shared_ptr<ofstream> excitonDistFile(new ofstream);
	excitonDistFile->open(excitonDistFileName);

	//Status updates
	int onePercent;
	if (autoComplete){ onePercent = 100; }
	else { onePercent = static_cast<int>(numSteps / 100.0); }
	int printCnt = 0;

	//For automatic simulation end there will be a comparison of differences to some threshold
	// If the differences, divided by the maximum current differences are < threshold for more
	// the numToFinish then the simulation will end.
	int numInARow = 0; //number of quotients that are below threshold in a row
	double difference = 0;//difference between average of numToCheck time step average num exciton points
	double maxDiff = 0; //The maximum difference between points used for difference.
	double prevAve =  numExcitonsAtCont; //average for last numToCheck time steps. Starts at init num of excitons
	double currAve = 0; //Average for current numToCheck time steps
	int timeSteps = 0; //current count of number of time steps
	bool simDone = false; //boolean symbolizing a finished simulation


	/*
	This section will consist of iterating until the maximum time has been reached. Each iteration
	for T will contain a single/multiple step for each exciton. 
	*/
	omp_set_num_threads(NUM_THREADS);
	while (T <= Tmax && !simDone) //iterates over time
	{
		//input loop
		while (_kbhit()) 
		{
			switch (_getch())
			{
				case 'Q':
					simDone = true;
			}
		}
		//reset the exciton count after each time step
		currCount = make_shared<vector<int>>(vector<int>(numRegions));
		T += deltaT; //set new time checkpoint
		#pragma omp parallel default(none) shared(CNT_List,currCount,regionBdr,gamma,excitons,deltaT)
		{
			#pragma omp for
			for (int exNum = 0; exNum < excitons->size(); exNum++) //iterates over excitons once
			{
				shared_ptr<exciton> currEx = (*excitons)[exNum];
				/*There is a change that the previous tr that was calculated was so long that it
				not only has extra time in the next deltaT, but it skips it completely. In this case
				the exciton movement is skipped all together and its extra time is decreased by deltaT
				until the exciton can move again. Otherwise the code runs as usual.
				*/
				double extraT = currEx->getTExtra();
				if (extraT > deltaT)
				{
					currEx->setTExtra(currEx->getTExtra() - deltaT);
					markCurrentExcitonPosition(CNT_List, currEx, currCount, regionBdr);
				}
				else
				{
					double tr_tot = extraT; //the sum of all tr's in the current deltaT time step
					/*give time to excitons that have no extra time. This happens only when excitons
					are just injected.*/
					if (extraT == 0)
					{
						tr_tot += -(1 / gamma)*log(getRand(true));
					}
					while (tr_tot <= deltaT)
					{
						//choose new state
						assignNextState(CNT_List, currEx, gamma, regionBdr);
						tr_tot += -(1 / gamma)*log(getRand(true)); // add the tr calculation to current time for individual particle
					}
					//Recording distribution
					markCurrentExcitonPosition(CNT_List, currEx, currCount, regionBdr);
					//record the time past deltaT that the assign next state will cover
					currEx->setTExtra(tr_tot - deltaT);
				}
			}
		}

		if (autoComplete)
		{
			timeSteps++;
			currAve += excitons->size();
			if (numToCheck == timeSteps)
			{
				currAve /= numToCheck; //calculate average
				//check for max difference
				if ((difference = currAve - prevAve) > maxDiff){ maxDiff = difference; }
				if (difference / maxDiff < threshold)
				{
					if (++numInARow >= numToFinish){ simDone = true; }
				}
				else{ numInARow = 0; }
				prevAve = currAve;
				timeSteps = 0; //reset the time counter
			}
		}

		//Output count vector to file since we want results after each time step.
		writeStateToFile(excitonDistFile, currCount, T);

		//Update Exciton List for injection and exit contact
		updateExcitonList(numExcitonsAtCont, excitons, currCount, inContact);
		//update time
		end = clock();
		runtime = diffclock(end, start);
		printCnt++;
		if (printCnt == onePercent)
		{
			ClearScreen();
			status = getRunStatus(T, Tmax, runtime, !autoComplete);
			cout << status << endl;
			printCnt = 0;
		}
		
	}

	//Close files finish program
	excitonDistFile->close();

	UINT64 numTSteps = static_cast<UINT64>(T / deltaT);

	//Write helper information
	writeExcitonDistSupportingInfo(tableFromFile, outputPath, numExcitonsAtCont, Tmax, deltaT, segLenMin, numRegions, xdim,
		minBin, rmax, numBins, lowAng, highAng, numAng, numTSteps, regLenMin, getRunTime(runtime));

	return 0;
}

/**
Removes excitons from the list if they are in the exit contact and inject excitons into the inContact if there are
not enough excitons

@param numExcitonAtCont The number of excitons at the injection contact
@param excitons The list of excitons to modify
@param currCount The count of excitons in each of the regions
@param inContact The list of segments in the injection contact
*/
void updateExcitonList(int numExcitonsAtCont, shared_ptr<vector<shared_ptr<exciton>>> excitons, 
	shared_ptr<vector<int>> currCount, shared_ptr<vector<shared_ptr<segment>>> inContact)
{
	// Adding excitons to the injection contact
	int numExAdd = 0; //The number of excitons to be added to the injection contact
	if ((numExAdd = numExcitonsAtCont - (*currCount)[0]) > 0)
	{
		for (int i = 0; i < numExAdd; i++)
		{
			excitons->push_back(make_shared<exciton>(exciton())); //initialize exciton at end of exciton list
			injectExciton((*excitons)[excitons->size()-1], inContact); //last element is the exciton to inject
		}
	}

	//Removing excitons from the exit contact
	int numExRem = (*currCount)[currCount->size() - 1]; //Last element in the count list = num of excitons to remove
	if (numExRem > 0)
	{
		for (int i = 0; i < excitons->size(); i++)
		{
			if ((*excitons)[i]->isAtOutContact())
			{
				shared_ptr<exciton> swap = (*excitons)[excitons->size() - 1];
				(*excitons)[excitons->size() - 1] = (*excitons)[i];
				(*excitons)[i] = swap;
				excitons->pop_back();
				//excitons->erase(excitons->begin() + i);
				i--;
			}
		}
	}
}


/**
Writes all passed information to a file that can be read by matlab to get the appropriate variable declarations
*/
void writeExcitonDistSupportingInfo(bool tableFromFile, string outputPath, int numExcitons, double Tmax, double deltaT, double segLenMin, int numRegions, 
	double xdim, double minBin, double rmax, int numBins, double lowAng, double highAng, int numAng, UINT64 numTSteps, double regLenMin, string runtime)
{
	string detailsFileName;
	if (tableFromFile)
	{
		detailsFileName = outputPath + "details_t.csv";
	}
	else
	{
		detailsFileName = outputPath + "details.csv";
	}
	shared_ptr<ofstream> detailsFile(new ofstream);
	detailsFile->open(detailsFileName);
	*detailsFile << "numExcitons,Tmax,deltaT,segLenMin,numRegions,xdim,minBin,rmax,numBins,lowAng,highAng,numAng,numTSteps,regLenMin,runtime" << endl;
	*detailsFile << numExcitons << "," << Tmax << "," << deltaT << "," << segLenMin << "," << numRegions << "," << xdim 
		<< "," << minBin << "," << rmax << "," << numBins << "," << lowAng << "," << highAng << "," << numAng << 
		"," << numTSteps << ","  << regLenMin << "," << runtime << endl;
}

/**
Writes the pertinent information of the current state of simulation to file

@param file The file to write to
@param currCount The distribution information to write to file
@param T The current time of the simulation
*/
void writeStateToFile(shared_ptr<ofstream> file, shared_ptr<vector<int>> currCount, double T)
{
	*file << T << "," << (*currCount)[0]; //Time followed by the counts
	for (UINT32 i = 1; i < currCount->size(); i++)
	{
		*file << "," << (*currCount)[i];
	}
	*file << endl;
}


/**
Takes the current exciton and adds it to the count of excitons in a certain region. 
This is called once per exciton per time step.

@param CNT_List The indexable list that contains the position information of exciton
@param exciton The exciton that we want to mark the position of.
@param currCount The vector that stores the count of excitons in the regions
@param regionBdr The region vector that allows selection of currCount index
*/
void markCurrentExcitonPosition(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> exciton, 
	shared_ptr<vector<int>> currCount, shared_ptr<vector<double>> regionBdr)
{
	//add count to the currCount vector in the location corresponding to the exciton 
	//   region location
	#pragma omp critical
	(*currCount)[getIndex(regionBdr,
		((*(*CNT_List)[exciton->getCNTidx()].segs)[exciton->getSegidx()]->mid(0)))]++;
}

/**
Checks to see if an exciton has moved into the out contact

@param exciton The exciton to check
@param regionBdr The array that stores contact information
@param CNT_List The indexable list that will reveal exciton location
@return True if the exciton has moved into the output contact, false otherwise
*/
bool hasMovedToOutContact(shared_ptr<exciton> exciton, shared_ptr<vector<double>> regionBdr, 
	shared_ptr<vector<CNT>> CNT_List)
{
	//If x component is creater than the second to last index of regionBdr, then it is in the output
	return (((*(*CNT_List)[exciton->getCNTidx()].segs)[exciton->getSegidx()]->mid(0)) 
					>= ((*regionBdr)[regionBdr->size() - 2]));
}


/**
Places the specified exciton into the input contact

@param exciton The exciton to add to the contact
@param inContact The list of segments for the injection contact
*/
void injectExciton(shared_ptr<exciton> exciton, shared_ptr<vector<shared_ptr<segment>>> inContact)
{
	energy E;
	if (static_cast<int>(round(getRand(false)))){ E = E11; }
	else{ E = E22; }
	exciton->setEnergy(E); //randomly set the energy of the exciton
	//choose a destination segment
	shared_ptr<segment> injectedSeg = (*inContact)[static_cast<int>(rand() % inContact->size())];
	//The self scattering table will have the correct indices for the current segment
	shared_ptr<vector<tableElem>> tbl = injectedSeg->tbl;
	//Set the exciton indices to the current segment
	exciton->setCNTidx((*tbl)[tbl->size() - 1].getTubeidx());
	exciton->setSegidx((*tbl)[tbl->size() - 1].getSegidx());
	exciton->setAtOutContact(false); //initialized out contact boolean

}


/**
Assigns the specified exciton to the next state in the simulation.

@param CNT_List The list of carbon nanotubes
@param e The exciton to be updated
@param inContact The list of segments for the injection contact
@param regionBdr The array that will help decide whether or not the exciton is in the out contact
*/
void assignNextState(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> e, double gamma, shared_ptr<vector<double>> regionBdr)
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
	e->setAtOutContact(hasMovedToOutContact(e, regionBdr, CNT_List));

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

@param vec The vector to search
@param prob The number to compare the indicies to
@return the index that requires the conditions in the method description
*/
int getIndex(shared_ptr<vector<double>> vec, double val)
{
	if (vec == nullptr)
	{
		ClearScreen();
		cout << "Error: Empty vector passed to getIndex()";
		system("pause");
		exit(EXIT_FAILURE);
	}
	int left = 0;
	int right = static_cast<int>(vec->size() - 1);
	while (left <= right)
	{
		if (right == left)
		{
			if ((*vec)[right] < val) { return right + 1; }
			return right;
		}
		int mid = static_cast<int>(static_cast<double>(right + left) / 2.0);
		double midVal = (*vec)[mid];
		if (midVal == val){ return mid; }

		if (midVal > val){ right = mid - 1; }
		else { left = mid + 1; }
	}
	return left;
}


/**
Takes a segment and determines which elements should be added to its tables based on distance away
from the segment.

@param CNT_List The list of carbon nanotubes to iterate over
@param seg The segment that we will be adding table elements to
@param maxDist If segments are within maxDist of seg, then they will be added to the table
@param colorMap A count of all rs at particular thetas to get mesh statistics
@param rs The r values that are needed to place values in the heat map
@param thetas The angles that are needed to place values in the heat map
@param addDataToTable Lambda to change how data is added to table based on build or read table.
@return The sum of all the rates calculated for the segment. For transition purposes
*/
double updateSegTable(double maxDist, dat2tab addDataToTable, tableUpdater &t, heatMapInfo &h)
{
	t.rate_tot = 0;
	//iterate over CNTs
	int i = 0; //CNT index counter
	//originally structured without i and j
	for (auto cntit = t.CNT_List->begin(); cntit != t.CNT_List->end(); ++cntit)
	{
		int j = 0; //segment index counter
		//iterate over all segments considered for seg
		for (auto segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			t.r = tableElem::calcDist((*t.seg)->mid, (*segit)->mid);
			t.theta = tableElem::calcThet(t.seg, segit);

			//Heat Map Additions
			if (t.r != 0)
			{
				(*h.map)[getIndex(h.rs, t.r)][getIndex(h.thetas, t.theta)]++; //////// Heat Map ///////
				//Check if within range
				if (t.r <= maxDist) /////// Building TABLE /////
				{
					t.i = i;
					t.j = j;
					addDataToTable(t);
				}
			}	
			j++;
		}
		i++;
	}
	return t.rate_tot;
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

/**
Measures time difference between two clocks

@param end end time
@param start start time
@return The difference in milliseconds
*/
double diffclock(clock_t end, clock_t start)
{

	double diffticks = end - start;
	double diffms = diffticks / (CLOCKS_PER_SEC / 1000);

	return diffms;
}

/**
Prints the status of the simulation

@param T The current time of the simulation
@param Tmax The end time of the simulation
@param runtime The amount of time the simulation has been running
@return The status of the simulation
*/
string getRunStatus(double T,double Tmax, double runtime, boolean runtimeKnown)
{
	string ret;
	if (runtimeKnown)
	{
		if (T != 0 && Tmax != 0)
		{
			ret = "Percent Complete: " + to_string(static_cast<int>(T / Tmax * 100))
				+ "%\n";
		}
		else
		{
			ret = "Preparing Simulation.\n";
		}
	}
	ret += "Time Simulated: " + to_string(T*1.0e9) + " ns\n";
	ret += "Runtime: ";
	ret += getRunTime(runtime);
	return ret;
}


/**
Gets the runtime of the simulation

@param runtime The total runtime of the simulation in milliseconds
@return A formated string to beter display run time
*/
string getRunTime(double runtime)
{
	string ret;

	int days = static_cast<int>(runtime / 86400000);
	if (days < 10){ ret += "0" + to_string(days) + ":"; }
	else{ ret += to_string(days) + ":"; }

	int hours = static_cast<int>(runtime / 3600000.) % 24;
	if (hours < 10){ ret += "0" + to_string(hours) + ":"; }
	else{ ret += to_string(hours) + ":"; }

	int minutes = static_cast<int>(runtime / 60000.) % 24;
	if (minutes < 10){ ret += "0" + to_string(minutes) + ":"; }
	else{ ret += to_string(minutes) + ":"; }

	int seconds = static_cast<int>(runtime / 1000.) % 60;
	if (seconds < 10){ ret += "0" + to_string(seconds); }
	else{ ret += to_string(seconds); }

	return ret;
}


/**
Clears the console
*/
void ClearScreen()
{
	HANDLE                     hStdOut;
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	DWORD                      count;
	DWORD                      cellCount;
	COORD                      homeCoords = { 0, 0 };

	hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
	if (hStdOut == INVALID_HANDLE_VALUE) return;

	/* Get the number of cells in the current buffer */
	if (!GetConsoleScreenBufferInfo(hStdOut, &csbi)) return;
	cellCount = csbi.dwSize.X *csbi.dwSize.Y;

	/* Fill the entire buffer with spaces */
	if (!FillConsoleOutputCharacter(
		hStdOut,
		static_cast<TCHAR>(' '),
		cellCount,
		homeCoords,
		&count
		)) return;

	/* Fill the entire buffer with the current colors and attributes */
	if (!FillConsoleOutputAttribute(
		hStdOut,
		csbi.wAttributes,
		cellCount,
		homeCoords,
		&count
		)) return;

	/* Move the cursor home */
	SetConsoleCursorPosition(hStdOut, homeCoords);
}

/**
Changes \ to / in strings. This is to fix file paths

@param path The path to be fixed
*/
string fixPath(string &path)
{
	regex rgx("\\\\");
	return regex_replace(path, rgx, "/");
}

//Returns the last Win32 error, in string format. Returns an empty string if there is no error.
string GetLastErrorAsString()
{
	//Get the error message, if any.
	DWORD errorMessageID = GetLastError();
	if (errorMessageID == 0)
		return string(); //No error message has been recorded

	LPSTR messageBuffer = nullptr;
	size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
		nullptr, errorMessageID, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), reinterpret_cast<LPSTR>(&messageBuffer), 0, nullptr);

	string message(messageBuffer, size);

	//Free the buffer.
	LocalFree(messageBuffer);

	return message;
}

/**
Function to add rates to table based on previous calculation

@param seg The segment the table element is being added to
@param r The separation between two segs
@param theta Angle between two segs
@param i CNT index
@param j segment index
@param rate_tot The running sum of the rates for the segment
@param c2c Unused. Only to fit function pointer requirements
*/
void addDataToTableCalc(tableUpdater &t)
{
	(*t.seg)->tbl->push_back(tableElem(t.r, t.theta, 6.4000e+19, t.i, t.j)); //tbl initialized in CNT::calculateSegments
	(*t.seg)->rateVec->push_back(t.rate_tot += ((*t.seg)->tbl->back()).getRate()); //tbl initialized in CNT::calculateSegments
}

/**
Function to add rates to table based on values read from a table

@param seg The segment the table element is being added to
@param r The separation between two segs
@param theta Angle between two segs
@param i CNT index
@param j segment index
@param rate_tot The running sum of the rates for the segment
@param c2c Amirhossein's tables used to extract rates
*/
void addDataToTableRead(tableUpdater &t)
{

}

/**
Uses selection sort algorithm to sort chiralities into an ordered list

@param list The list of Chiralities to sort
*/
void sortChiralities(vector<Chirality> &list)
{
	for (int i = 0; i < list.size()-1; i++)
	{
		int minIdx = i;
		Chirality minChir = list[minIdx];
		for (int j = i+1; j < list.size(); j++)
		{
			if (minChir.compare(list[j]) > 0)
			{
				minIdx = j;
				minChir = list[j];
			}
		}
		Chirality temp = list[i];
		list[i] = minChir;
		list[minIdx] = temp;
	}
}

/**
Given the chirality, finds the index at which that chirality belongs

@param vec The vector of chiralities to search through
@param val The chirality which index we wish to find
@return Index that the chirality belongs to. -1 if fail or error.
*/
int getIndex(vector<Chirality> &vec, Chirality &val)
{
	int left = 0;
	int right = static_cast<int>(vec.size() - 1);
	int center;

	while (left <= right)
	{
		if (right == left)
		{
			if (vec[right] == val){ return right; }
			return -1;
		}
		center = (right + left) / 2;
		if (vec[center] == val){ return center; }

		if (vec[center].compare(val) > 0){ right = center - 1; }
		else { left = center + 1; }
	}
	return -1;
}

/**
Adds chiralities to list from file name

@param vec The vector to add the chiralities to
@param filename The filename to parse chiralities from
*/
void addChiralitiesToList(vector<Chirality> &vec, char* filename)
{
	chirPair pair = getChiralityFromFilename(filename);
	if (find(vec.begin(), vec.end(), pair.c1) == vec.end())
	{
		vec.push_back(pair.c1);
	}
	if (find(vec.begin(), vec.end(), pair.c2) == vec.end())
	{
		vec.push_back(pair.c2);
	}
}


/**
Gets chirality information from the binary file names
*/
chirPair getChiralityFromFilename(char* filename)
{
	chirPair ret_pair;
	vector<int> intMatch = vector<int>(4);
	char *next_token = nullptr;
	char* grabToken = strtok_s(filename, ",_.", &next_token);
	for (int i = 0; i < 4; i++)
	{
		intMatch[i] = atoi(grabToken);
		grabToken = strtok_s(nullptr, ",_.", &next_token);
	}
	ret_pair.c1 = Chirality(intMatch[0], intMatch[1]);
	ret_pair.c2 = Chirality(intMatch[2], intMatch[3]);

	return ret_pair;
}
