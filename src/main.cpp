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
#include <omp.h>
#include <stdint.h>
#include <unistd.h> // used for changing the working directory in the program: chdir

#include "file_management.h"
#include "functions.h"
#include "write_log.h"


using namespace std;

#define MAX_T_STEPS 10000000 //maximum number of time steps allowed for the simulation

//Global variables
double ymax = 0; //stores maximum height the cylinders of the CNTs are found at. All will be greater than 0.

//Runs the file input, monte carlo, and file output sections of code
int main(int argc, char *argv[])
{

	// unsigned int NUM_THREADS = omp_get_max_threads();
	string status;
	//Varible initialization
	double segment_length = 100.0; //[Angstroms]
	double regLenMin = segment_length;
	double runtime;
	int numSteps = 10000; //number of deltaT's the simulation runs if the user specifies
	double maxDist = 300; //[Angstroms] The maximum segment separation to be considered for exciton transfer
	//The number of excitons to be in injection contact at start
	int numExcitonsAtCont = 0;
	bool autoComplete = false;
	double threshold = .01; //The value that the difference/maxDiff must be below to finish
	int numToFinish = 5;// Number of differences/maxDiff  that must be below thresh in a row to finish
	int numToCheck = 1000; //number of time steps to check finish completion
	double tfac = log(.3);

	//timing functions
	clock_t start = clock();

	//Initialize random number generator before anything to ensure that getRand() always works
	initRandomNumGen();

	file_management file_manager;

	if (argc == 2)
	{
		file_manager.input_directory = argv[1];
		// int rc = chdir(file_manager.input_directory.c_str());
		// if (rc < 0)
		// {
		// 	cout << "unable to change the output directory!!!" << endl;
		// 	exit(EXIT_FAILURE);
		// }
	}
	else
	{
		cout << "input directory must be entered as an argument!!!" << endl;;
		exit(EXIT_FAILURE);
	}


	// string log_input = "test!!!!";
	// cout << log_input << endl;
	// write_log(2);

	// return 0;

	////////////////////////////////// INPUT FOLDER ///////////////////////////////


	// check if folder or file exists.
	struct stat buf;
	if (stat(file_manager.input_directory.c_str(),&buf) != 0)
	{
		cout << "Error: incorrect result directory path!!!" << endl;;
		exit(EXIT_FAILURE);
	}

	{
		char last_char = file_manager.input_directory.at(file_manager.input_directory.size()-1);
		if (last_char != '/')
		{
			file_manager.input_directory.push_back('/');
		}
		cout << file_manager.input_directory << endl;
	}

	//////////////////////////// BUILD FILE LIST ///////////////////////////////////////////
	// This next block creates a list of files in the directory chosen to be looked at.

	DIR *resDir;
	struct dirent *ent = nullptr;
	unique_ptr<list<string>>  fileList(new list<string>(0));

	regex rgx("CNT_Num_\\d+\\.csv"); //files we want to look through

	//Check if folder can be opened - should work due to above checks
	if ((resDir = opendir(file_manager.input_directory.c_str())) != nullptr)
	{
		//throw away first two results as they are . and ..
		readdir(resDir);
		readdir(resDir);

		//iterate over all of the real files
		while ((ent = readdir(resDir)) != nullptr)
		{
			smatch matches; //match_results for string objects
			string tmps =  string(ent->d_name);
			regex_search(tmps, matches, rgx);
			if (!matches.empty())
			{
				fileList->push_back(ent->d_name);
			}
		}
		closedir(resDir); //deletes pointer

	}
	else
	{
		cout << "Could not open results directory!!!" << endl;
		exit(EXIT_FAILURE);
	}
	delete ent;
	ent = nullptr;

	/////////////////////////////////// XML FILE PARSE /////////////////////////////////////
	
	double rmax = 0; //maximum possible differences between sections of CNTs
	double xdim = 0; //Dimension in which exciton populations will be monitored
	double ydim = 0;
	double zdim = 0;

	{

		string inputXMLPath = file_manager.input_directory + "input.xml";

		rapidxml::xml_document<> doc; //create xml object
		rapidxml::file<> xmlFile(inputXMLPath.c_str()); //open file
		doc.parse<0>(xmlFile.data()); //parse contents of file
		rapidxml::xml_node<>* currNode = doc.first_node(); //gets the node "Document" or the root node


		currNode = currNode->first_node("outputDirectory");
		file_manager.output_directory = currNode->value();

		cout << "output_directory = " << file_manager.output_directory << endl;

		// DEVICE DIMENSIONS NODE //
		currNode = currNode->next_sibling("DeviceDimensions"); //Output folder

		xdim = convertUnits(string(currNode->first_node("units")->value()),	atof(currNode->first_node("xdim")->value()));
		ydim = convertUnits(string(currNode->first_node("units")->value()),	atof(currNode->first_node("ydim")->value()));
		zdim = convertUnits(string(currNode->first_node("units")->value()), atof(currNode->first_node("zdim")->value()));
		//incorrect units
		if (xdim == INT_MIN || ydim == INT_MIN || zdim == INT_MIN)
		{
			cout << "Error: Incorrect units for device dimensions!!!" << endl;
			exit(EXIT_FAILURE);
		}
		//incorrect range
		else if (xdim <= 0 || ydim <= 0 || zdim <= 0)
		{
			cout << "Error: Must enter positive device dimensions!!!" << endl;
			exit(EXIT_FAILURE);
		}

		cout << "xdim = " << xdim << endl;
		cout << "ydim = " << ydim << endl;
		cout << "zdim = " << zdim << endl;
		// END DEVICE DIMENSIONS NODE //

		//x and z dim are the ground dimensions and y is height. Need to make sure bottom corner
		// to top other corner is included in r range. ymax is not defined so this is intermediate
		// value. It is changed after CNT initialization.
		rmax = sqrt(pow(xdim,2)+pow(zdim,2));

		// NUMBER OF EXCITONS NODE //
		currNode = currNode->next_sibling("numberExcitons"); //to number of excitons
		numExcitonsAtCont = atoi(currNode->value());
		if (numExcitonsAtCont <= 0)
		{
			cout << "Error: Must have positive number of excitons." << endl;
			exit(EXIT_FAILURE);
		}
		cout << "numExcitonsAtCont = " << numExcitonsAtCont << endl;
		// END NUMBER OF EXCITONS NODE //


		// REGION LENGTH NODE //
		currNode = currNode->next_sibling("regionLength");
		regLenMin = convertUnits(string(currNode->first_node("units")->value()), atof(currNode->first_node("min")->value()));
		if (regLenMin <= 0 || regLenMin > xdim)
		{
			cout << "Error: Region length must be positive number!!!" << endl;
			exit(EXIT_FAILURE);
		}
		// END REGION LENGTH NODE //

		// SEGMENT LENGTH NODE //
		currNode = currNode->next_sibling("segmentLength");
		segment_length = convertUnits(string(currNode->first_node("units")->value()), atof(currNode->first_node("min")->value()));
		if (segment_length <= 0)
		{
			cout << "Error: Segment length must be positive number!!!" << endl;
			exit(EXIT_FAILURE);
		}

		cout << "segment_length = " << segment_length << endl;
		// END SEGMENT LENGTH NODE //

		// SEGMENT SEPARATION //
		currNode = currNode->next_sibling("segmentSeparation");
		maxDist = convertUnits(string(currNode->first_node("units")->value()), atof(currNode->first_node("min")->value()));
		if (maxDist <= 0)
		{
			cout << "Error: Segment separation must be a positive number!!!" << endl;
			exit(EXIT_FAILURE);
		}
		// END SEGMENT SEPARATION //

		// NUMBER OF TIME STEPS //
		currNode = currNode->next_sibling("numberTimeSteps");
		numSteps = atoi(currNode->value());
		if (numSteps <= 0)
		{
			cout << "Error: Must have positive number time steps.!!!" << endl;
			exit(EXIT_FAILURE);
		}
		// END NUMBER OF TIME STEPS //


		// PERCENT FREE FLIGHT TIMES ABOVE DELTA T //
		currNode = currNode->next_sibling("percentFreeFlightTimesAboveDeltaT");
		tfac = atoi(currNode->value()) / 100.0;
		if (tfac <= 1 && tfac > 0)
		{
			tfac = abs(log(tfac));
		}
		else
		{
			cout << "Error: Percent of free flight times above delta T must be greater than 0 and less than or equal to 100!!!" << endl;
			exit(EXIT_FAILURE);
		}
		// END PERCENT FREE FLIGHT TIMES ABOVE DELTA T  //

		// AUTO COMPLETE //
		currNode = currNode->next_sibling("autoComplete")->first_node("enabled");

		if (string(currNode->value()).compare("true") == 0)
		{
			autoComplete = true;
			numSteps = MAX_T_STEPS;

			currNode = currNode->next_sibling("threshold");
			threshold = atof(currNode->value());
			if (threshold < 0.0 || threshold > 1.0)
			{
				cout << "Configuration Error: Threshold must be value greater than 0 and less than 1!!!" << endl;
				exit(EXIT_FAILURE);
			}

			currNode = currNode->next_sibling("numBelowThreshold");
			numToFinish = atoi(currNode->value());
			if (numToFinish <= 0) //input validation
			{
				cout << "Configuration Error: Number below threshold must be greater than 0!!!" << endl;
				exit(EXIT_FAILURE);
			}

			int numToCheck = atoi(currNode->value());
			if (numToCheck <= 0)
			{
				cout << "Configuration Error: Number to average must be greater than 0!!!" << endl;
				exit(EXIT_FAILURE);
			}
		}
		// END AUTO COMPLETE // 
	}

	//@@@@@@@@ TIME UPDATE @@@@@@@@//
	clock_t end = clock();
	runtime = diffclock(end, start);
	status = getRunStatus(0, 0, runtime, !autoComplete);
	cout << status << endl;
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@//

	return 0;

	/*
	This section creates the list of CNTs with all of the relevant information provided in their
	respective files getting processed.
	*/

	//////////////////////////// BUILD CNT AND SEGMENTS //////////////////////////////////////////

	//Iterate through the files and extract
	shared_ptr<vector<CNT>> CNT_List(new vector<CNT>(0));

	for (list<string>::iterator it = fileList->begin(); it != fileList->end(); ++it)
	{

		CNT_List->push_back(CNT(*(it), file_manager.input_directory, segment_length));
	}

	//Extra check to ensure that all initilizations were successful
	for (vector<CNT>::iterator it = CNT_List->begin(); it != CNT_List->end(); ++it)
	{
		if (!(*it).isInitialized())
		{
			cout << "CNT initialization failure!!!" << endl;
			exit(EXIT_FAILURE);
		}
	}

	cout << "the cnts are read correctly!!!" << endl;

	//Set final rmax value
	rmax = sqrt(pow(rmax, 2) + pow(ymax,2));


	//@@@@@@@@ TIME UPDATE @@@@@@@@//
	end = clock();
	runtime = diffclock(end, start);
	// ClearScreen();
	status = getRunStatus(0, 0, runtime, !autoComplete);
	cout << status << endl;
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@//

	// return 0;

	//////////////////////////// CREATE DIST VECTORS AND ARRAYS ///////////////////////////////

	rmax = 100 * ceil(rmax / 100.0);
	int numBins = static_cast<int>(rmax / 10.0); //number of bins to place r's into 
	double minBin = rmax / static_cast<double>(numBins); //[Angstroms] The size of the bins
	shared_ptr<vector<double>> rs = linspace(minBin, rmax, numBins); //Builds rs vector within valid r range

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
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		double newGamma;
		//loop over segments in each CNTs
		for (vector<shared_ptr<segment>>::iterator segit = cntit->segments->begin(); segit != cntit->segments->end(); ++segit)
		{
			int regIdx = getIndex(regionBdr, (*segit)->mid(0));
			(*segCountPerReg)[regIdx]++; //increment the count based on where segment is
			if (regIdx == 0){ inContact->push_back(*segit); } //First region is injection contact
			//get add to each segment relevant table entries
			newGamma = updateSegTable(CNT_List, segit, maxDist, heatMap, rs, thetas);
			if (newGamma > gamma){ gamma = newGamma; }
			numSegs++;
		}
		buildTblCntr++;
		if (buildTblCntr == buildTblDecPerc)
		{
			//@@@@@@@@ TIME UPDATE @@@@@@@@//
			end = clock();
			runtime = diffclock(end, start);
			// ClearScreen();
			status = getRunStatus(0, 0, runtime, !autoComplete);
			cout << status << endl;
			//@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
			buildTblCntr = 0;
		}
		
	}


	/////////////////////////////// OUTPUT SEG COUNT PER REGION //////////////////////////////////////

	string segmentCountFileName = file_manager.output_directory + "segmentCountPerRegion.csv";
	ofstream segmentCountFile;
	segmentCountFile.open(segmentCountFileName);
	segmentCountFile << (*segCountPerReg)[0];
	//iterate through all of the rs and then thetas while printing to file
	for (uint32_t i = 1; i < segCountPerReg->size(); i++)
	{
		segmentCountFile << "," << (*segCountPerReg)[i];
	}
	segmentCountFile.close();

	/////////////////////////////// OUTPUT HEATMAP //////////////////////////////////////

	string heatMapFileName = file_manager.output_directory + "heatMap.csv";
	ofstream heatMapFile;
	heatMapFile.open(heatMapFileName);
	//iterate through all of the rs and then thetas while printing to file
	for (uint32_t i = 0; i < rs->size(); i++)
	{
		for (uint32_t j = 0; j < thetas->size(); j++)
		{
			heatMapFile << (*heatMap)[i][j] << ",";
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
	string excitonDistFileName = file_manager.output_directory + "excitonDist.csv";
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
	// omp_set_num_threads(NUM_THREADS);
	while (T <= Tmax && !simDone) //iterates over time
	{

		// break the simulation if the input key is "Q"
		//input loop
		// while (_kbhit()) 
		// {
		// 	switch (_getch())
		// 	{
		// 		case 'Q':
		// 			simDone = true;
		// 	}
		// }

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
			// ClearScreen();
			status = getRunStatus(T, Tmax, runtime, !autoComplete);
			cout << status << endl;
			printCnt = 0;
		}
		
	}

	//Close files finish program
	excitonDistFile->close();

	uint64_t numTSteps = static_cast<uint64_t>(T / deltaT);

	//Write helper information
	writeExcitonDistSupportingInfo(file_manager.output_directory, numExcitonsAtCont, Tmax, deltaT, segment_length, numRegions, xdim,
		minBin, rmax, numBins, lowAng, highAng, numAng, numTSteps, regLenMin, getRunTime(runtime));

	return 0;
}

