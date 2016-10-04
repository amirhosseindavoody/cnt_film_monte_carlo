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

#include "functions.h"


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
void writeExcitonDistSupportingInfo(string outputPath, int numExcitons, double Tmax, double deltaT, double segLenMin, int numRegions, 
	double xdim, double minBin, double rmax, int numBins, double lowAng, double highAng, int numAng, uint64_t numTSteps, double regLenMin, string runtime)
{
	string detailsFileName = outputPath + "details.csv";
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
	for (uint32_t i = 1; i < currCount->size(); i++)
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
	exciton->setEnergy(static_cast<int>(round(getRand(false)) + 1)); //randomly set the energy of the exciton
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

			//Heat Map Additions
			if (r != 0)
			{
				(*heatMap)[getIndex(rs, r)][getIndex(thetas, theta)]++; //////// Heat Map ///////
				//Check if within range
				if (r <= maxDist) /////// Building TABLE /////
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
string getRunStatus(double T,double Tmax, double runtime, bool runtimeKnown)
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