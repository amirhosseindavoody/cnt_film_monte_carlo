// CNTMC_Project.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <list>
#include <vector>
#include <memory>
#include <locale>
#include <math.h>
#include <regex>
#include <omp.h>
#include <stdint.h>
#include <unistd.h> // used for changing the working directory in the program: chdir

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "file_management.h"
#include "functions.h"
#include "write_log.h"
#include "CNT.h"
#include "tableElem.h"
#include "simulation_parameters.h"


using namespace std;

//Global variables
double ymax = 0; //stores maximum height the cylinders of the CNTs are found at. All will be greater than 0.

//Runs the file input, monte carlo, and file output sections of code
int main(int argc, char *argv[])
{
	double runtime;
	clock_t start = clock(); //timing functions

	//Initialize random number generator before anything to ensure that getRand() always works
	init_random_number_generator();

	if (argc != 2)
	{
		cout << "input directory must be entered as an argument!!!" << endl;;
		exit(EXIT_FAILURE);
	}

	file_management file_manager;

	file_manager.set_input_directory(argv[1]);
	file_manager.build_cnt_file_list();

	simulation_parameters sim = file_manager.parse_xml();

	file_manager.change_working_directory(file_manager.output_directory.c_str());

	//////////////////////////// BUILD CNT AND SEGMENTS //////////////////////////////////////////

	//Iterate through the files and extract
	shared_ptr<vector<CNT>> CNT_List = make_shared<vector<CNT>>();

	for (list<string>::iterator curr_file = file_manager.file_list->begin(); curr_file != file_manager.file_list->end(); ++curr_file)
	{
		string file_name = *curr_file;
		string path = file_manager.get_input_directory();
		double segment_length = sim.segment_length;

		CNT_List->push_back(CNT(file_name, path, segment_length));
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

	//Set final sim.rmax value
	sim.rmax = sqrt(pow(sim.rmax, 2) + pow(ymax,2));

	// return 0;

	//////////////////////////// CREATE DIST VECTORS AND ARRAYS ///////////////////////////////

	sim.rmax = 100 * ceil(sim.rmax / 100.0);
	int numBins = static_cast<int>(sim.rmax / 10.0); //number of bins to place r's into 
	double minBin = sim.rmax / static_cast<double>(numBins); //[Angstroms] The size of the bins
	shared_ptr<vector<double>> rs = linspace(minBin, sim.rmax, numBins); //Builds rs vector within valid r range

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
	int numRegions = static_cast<int>(sim.xdim / sim.region_length_min); //number of regions in the simulation
	double extra = sim.xdim - numRegions*sim.region_length_min; //extra amount * numRegions that should be added to sim.region_length_min
	double regLen = sim.region_length_min + extra / numRegions; //The length of each region.

	//The boundary of the rgions in the x direction
	shared_ptr<vector<double>> regionBdr = linspace(regLen - (sim.xdim / 2), sim.xdim / 2, numRegions);
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
			newGamma = updateSegTable(CNT_List, segit, sim.maximum_distance, heatMap, rs, thetas);
			if (newGamma > gamma){ gamma = newGamma; }
			numSegs++;
		}
		buildTblCntr++;
		if (buildTblCntr == buildTblDecPerc)
		{
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
	shared_ptr<vector<shared_ptr<exciton>>> excitons(new vector<shared_ptr<exciton>>(sim.number_of_excitons));
	for (int exNum = 0; exNum < sim.number_of_excitons; exNum++)
	{
		(*excitons)[exNum] = make_shared<exciton>(exciton()); //initialize exciton at location in exciton list
		injectExciton((*excitons)[exNum], inContact);
	}


	//////////////////////////////////// TIME STEPS ///////////////////////////////////
	double deltaT = (1 / gamma)*sim.tfac; //time steps at which statistics are calculated
	double Tmax = deltaT * sim.number_of_steps; //maximum simulation time
	double T = 0; //Current simulation time, also time at which next stats will be calculated

	shared_ptr<vector<int>> currCount; //The count of excitons in each region of the CNT mesh

	//File output initializations
	string excitonDistFileName = file_manager.output_directory + "excitonDist.csv";
	shared_ptr<ofstream> excitonDistFile(new ofstream);
	excitonDistFile->open(excitonDistFileName);

	//Status updates
	int onePercent;
	if (sim.auto_complete){ onePercent = 100; }
	else { onePercent = static_cast<int>(sim.number_of_steps / 100.0); }
	int printCnt = 0;

	//For automatic simulation end there will be a comparison of differences to some sim.threshold
	// If the differences, divided by the maximum current differences are < sim.threshold for more
	// the sim.num_to_finish then the simulation will end.
	int numInARow = 0; //number of quotients that are below sim.threshold in a row
	double difference = 0;//difference between average of sim.num_to_check time step average num exciton points
	double maxDiff = 0; //The maximum difference between points used for difference.
	double prevAve =  sim.number_of_excitons; //average for last sim.num_to_check time steps. Starts at init num of excitons
	double currAve = 0; //Average for current sim.num_to_check time steps
	int timeSteps = 0; //current count of number of time steps
	bool simDone = false; //boolean symbolizing a finished simulation


	/*
	This section will consist of iterating until the maximum time has been reached. Each iteration
	for T will contain a single/multiple step for each exciton. 
	*/
	// omp_set_num_threads(NUM_THREADS);
	while (T <= Tmax && !simDone) //iterates over time
	{

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

		if (sim.auto_complete)
		{
			timeSteps++;
			currAve += excitons->size();
			if (sim.num_to_check == timeSteps)
			{
				currAve /= sim.num_to_check; //calculate average
				//check for max difference
				if ((difference = currAve - prevAve) > maxDiff){ maxDiff = difference; }
				if (difference / maxDiff < sim.threshold)
				{
					if (++numInARow >= sim.num_to_finish){ simDone = true; }
				}
				else{ numInARow = 0; }
				prevAve = currAve;
				timeSteps = 0; //reset the time counter
			}
		}

		//Output count vector to file since we want results after each time step.
		writeStateToFile(excitonDistFile, currCount, T);

		//Update Exciton List for injection and exit contact
		updateExcitonList(sim.number_of_excitons, excitons, currCount, inContact);

		printCnt++;
		if (printCnt == onePercent)
		{
			printCnt = 0;
		}
		
	}

	//Close files finish program
	excitonDistFile->close();

	// uint64_t numTSteps = static_cast<uint64_t>(T / deltaT);

	//Write helper information
	// writeExcitonDistSupportingInfo(file_manager.output_directory, sim.number_of_excitons, Tmax, deltaT, sim.segment_length, numRegions, sim.xdim, minBin, sim.rmax, numBins, lowAng, highAng, numAng, numTSteps, sim.region_length_min, getRunTime(runtime));

	return 0;
}

