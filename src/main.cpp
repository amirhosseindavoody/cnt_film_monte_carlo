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
	vector<CNT> cnt_list;
	for (int i=0; i < file_manager.file_list.size(); i++)
	{
		string file_name = file_manager.file_list[i];
		string path = file_manager.get_input_directory();
		double segment_length = sim.segment_length;

		cnt_list.push_back(CNT(file_name, path, segment_length));
	}

	//////////////////////////// CREATE DIST VECTORS AND ARRAYS ///////////////////////////////

	sim.rmax = sqrt(pow(sim.rmax, 2) + pow(ymax,2));

	int n_r = 100;
	int n_theta = 100;
	vector<double> rs = linspace(0, sim.rmax, n_r);
	vector<double> thetas = linspace(0, 90, n_theta);
	vector<vector<int>> heatMap(n_r, vector<int>(n_theta));
	

	//////////////////////////// BUILD TABLE ////////////////////////////////////////////////

	///////////// Variables for placing excitons in the future /////////////
	int numRegions = static_cast<int>(sim.xdim / sim.region_length_min); //number of regions in the simulation
	double extra = sim.xdim - numRegions*sim.region_length_min; //extra amount * numRegions that should be added to sim.region_length_min
	double regLen = sim.region_length_min + extra / numRegions; //The length of each region.

	vector<double> regionBdr = linspace(regLen - (sim.xdim / 2), sim.xdim / 2, numRegions); //The boundary of the rgions in the x direction
	vector<int> segCountPerReg(numRegions); //To get dist stats
	vector<shared_ptr<segment>> inContact(0); //List of segments in the first region, which is used as a input contact

	
	//iterate through all of the CNTs and segments
	int numSegs = 0; //The total number of segments in the simulation, used in exciton placement
	double gamma = 0; //The maximum of sums of gammas from each segment. This sets the constant gamma value for the entire simulation

	int buildTblDecPerc = static_cast<int>(cnt_list.size() / 10.0); //ten percent of total number of tubes
	int buildTblCntr = 0;
	//loop over CNTs
	for (int i=0; i<cnt_list.size(); i++)
	{
		CNT &curr_cnt = cnt_list[i];

		double newGamma;
		
		//loop over segments in each CNTs
		for (int j=0; j<curr_cnt.segments.size(); j++)
		{
			segment &curr_segment = curr_cnt.segments[j];

			int regIdx = getIndex(regionBdr, curr_segment.mid(0));
			segCountPerReg[regIdx]++; //increment the count based on where segment is

			if (regIdx == 0){ inContact.push_back(make_shared<segment>(curr_segment)); } //First region is injection contact
			//get add to each segment relevant table entries
			newGamma = updateSegTable(cnt_list, curr_segment, sim.maximum_distance, heatMap, rs, thetas);
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
	for (uint32_t i = 0; i < segCountPerReg.size(); i++)
	{
		segmentCountFile << " , " << segCountPerReg[i];
	}
	segmentCountFile.close();

	/////////////////////////////// OUTPUT HEATMAP //////////////////////////////////////

	string heatMapFileName = file_manager.output_directory + "heatMap.csv";
	ofstream heatMapFile;
	heatMapFile.open(heatMapFileName);
	//iterate through all of the rs and then thetas while printing to file
	for (uint32_t i = 0; i < rs.size(); i++)
	{
		for (uint32_t j = 0; j < thetas.size(); j++)
		{
			heatMapFile << heatMap[i][j] << " , ";
		}
		heatMapFile << "\n";
	}
	heatMapFile.close();



	/////////////////////////////// ADD SELF SCATTERING //////////////////////////////////

	addSelfScattering(cnt_list, gamma);

	/*
	Now that the tables have been built, the next step is to populate the mesh with excitons. The way this will happen
	is, after a specific number of excitons are chosen to be created, each exciton will be created and assigned an index
	that corresponds to the CNT mesh. The available locations for excitons to start will be determined by the number of 
	segments that are within the delta x region specified. The exciton assignment will be equally and randomly distributed 
	among the list of segments in the delta x region.
	*/

	//////////////////////////// PLACE EXCITONS ////////////////////////////////////////////

	vector<exciton> excitons(sim.number_of_excitons);
	for (int exNum = 0; exNum < sim.number_of_excitons; exNum++)
	{
		injectExciton(excitons[exNum], inContact);
	}


	//////////////////////////////////// TIME STEPS ///////////////////////////////////
	double deltaT = (1 / gamma)*sim.tfac; //time steps at which statistics are calculated
	double Tmax = deltaT * sim.number_of_steps; //maximum simulation time
	double T = 0; //Current simulation time, also time at which next stats will be calculated

	vector<int> currCount(numRegions); //The count of excitons in each region of the CNT mesh

	//File output initializations
	string excitonDistFileName = file_manager.output_directory + "excitonDist.csv";
	ofstream excitonDistFile;
	excitonDistFile.open(excitonDistFileName);

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

	// monte carlo stepping loop
	while (T <= Tmax && !simDone) //iterates over time
	{
		T += deltaT;

		//reset the exciton count after each time step
		clear_vector(currCount, 0);

		for (int exNum = 0; exNum < excitons.size(); exNum++)
		{
			exciton &curr_exciton = excitons[exNum];
			/*There is a change that the previous tr that was calculated was so long that it
			not only has extra time in the next deltaT, but it skips it completely. In this case
			the exciton movement is skipped all together and its extra time is decreased by deltaT
			until the exciton can move again. Otherwise the code runs as usual.
			*/
			double extraT = curr_exciton.getTExtra();
			if (extraT > deltaT)
			{
				curr_exciton.setTExtra(extraT - deltaT);
				markCurrentExcitonPosition(cnt_list, curr_exciton, currCount, regionBdr);
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
					assignNextState(cnt_list, curr_exciton, gamma, regionBdr);
					tr_tot += -(1 / gamma)*log(getRand(true)); // add the tr calculation to current time for individual particle
				}
				//Recording distribution
				markCurrentExcitonPosition(cnt_list, curr_exciton, currCount, regionBdr);
				//record the time past deltaT that the assign next state will cover
				curr_exciton.setTExtra(tr_tot - deltaT);
			}
		}

		if (sim.auto_complete)
		{
			timeSteps++;
			currAve += excitons.size();
			if (sim.num_to_check == timeSteps)
			{
				currAve /= sim.num_to_check; //calculate average
				
				//check for max difference
				if ((difference = currAve - prevAve) > maxDiff)
				{
					maxDiff = difference;
				}
				
				if (difference / maxDiff < sim.threshold)
				{
					if (++numInARow >= sim.num_to_finish)
					{
						simDone = true;
					}
				}
				else
				{
					numInARow = 0;
				}
				prevAve = currAve;
				timeSteps = 0; //reset the time counter
			}
		}

		//Output count vector to file since we want results after each time step.
		writeStateToFile(excitonDistFile, currCount, T);

		//Update Exciton List for injection and exit contact
		updateExcitonList(sim.number_of_excitons, excitons, currCount, inContact);
		
	}

	//Close files finish program
	excitonDistFile.close();

	return 0;
}