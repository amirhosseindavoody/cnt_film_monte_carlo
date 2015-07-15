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


using namespace std;

#define TFAC 2.302585092994046 //Factor to multiply 1/gamma by to get 90% of rand numbers giving tr less than deltaT

//method declarations
string folderPathPrompt(bool incorrect);
double updateSegTable(shared_ptr<vector<CNT>> CNT_List, vector<segment>::iterator seg, double maxDist);
int getIndex(shared_ptr<vector<double>> vec, double prob);
int getIndex(shared_ptr<vector<double>> vec, double prob, int left, int right);
double getRand();
void addSelfScattering(shared_ptr<vector<CNT>> CNT_List, double maxGam);
void assignNextState(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> e, double gamma);

int main(int argc, char *argv[])
{

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

	//Checking validity of folder path
	bool done = false;
	while (!done)
	{
		//Container for file/folder information
		struct stat info;

		if (stat(resultFolderPath.c_str(), &info) != 0)
		{
			printf("Cannot access %s\n", resultFolderPath.c_str());
			resultFolderPath = folderPathPrompt(true);
		}
		else if (info.st_mode & S_IFDIR)
		{
			done = true;
		}
		else
		{
			printf("%s is not a directory\n", resultFolderPath.c_str());
			resultFolderPath = folderPathPrompt(true);
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
		readdir(resDir);readdir(resDir);
		//iterate over all of the real files
		while ((ent = readdir(resDir)) != nullptr)
		{
			fileList->push_back(ent->d_name);
		}
		closedir(resDir); //deletes pointer
	} else
	{
		cout << "Could not open directory. Please try program again.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	delete ent;


	/*
	This section creates the list of CNTs with all of the relevant information provided in their
	respective files getting processed.
	*/

	//////////////////////////// BUILD CNT AND SEGS //////////////////////////////////////////

	//Iterate through the files and extract
	shared_ptr<vector<CNT>> CNT_List(new vector<CNT>(0));

	for (list<string>::iterator it = fileList->begin(); it != fileList->end(); ++it)
	{
		CNT_List->push_back(CNT(*(it), resultFolderPath,100.0));
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
		for (vector<segment>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			//get add to each segment relelant table entries
			newGamma = updateSegTable(CNT_List, segit, maxDist);
			if (newGamma > gamma){ gamma = newGamma; }
			numSegs++;
		}
	}

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

	//Initialize random number generation
	time_t seconds;
	time(&seconds); //assign time from clock
	//Seed the random number generator
	srand(static_cast<int>(seconds));

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
			(*excitons)[exNum]->setEnergy(static_cast<int>(round(getRand()) + 1));
			(*excitons)[exNum]->setCNTidx(getIndex(cntProb, getRand()));
			//sets seg index to a number between 0 and the number of segs for the cnt - 1.
			//  The complicated part gets the number of segments from the current cnt
			(*excitons)[exNum]->setSegidx(rand() % ((*CNT_List)[(*excitons)[exNum]->getCNTidx()].segs->size()));
			shared_ptr<vector<segment>> temp((*CNT_List)[(*excitons)[exNum]->getCNTidx()].segs); //get segment list
			done = (*temp)[(*excitons)[exNum]->getSegidx()].setExciton((*excitons)[exNum]);
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
				tr_tot+= -(1 / gamma)*log(getRand()); // add the tr calculation to current time for individual particle
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
	bool done = false; //flag for completion
	//Segment the current exciton is located on
	segment seg = (*((*CNT_List)[e->getCNTidx()].segs))[e->getSegidx()]; 
	while (!done)
	{
		int tblIdx = getIndex(seg.rateVec, getRand()*gamma);
		tableElem tbl = (*seg.tbl)[tblIdx];
		segment newSeg = (*((*CNT_List)[tbl.getTubeidx()].segs))[tbl.getSegidx()];
		done = newSeg.setExciton(e);
		if (done)
		{
			if (!seg.removeExciton(e))
			{
				cout << "Error: Exciton could not be removed from segment.\n";
				system("pause");
				exit(EXIT_SUCCESS);
			}
			e->setCNTidx(tbl.getTubeidx());
			e->setSegidx(tbl.getSegidx());
		}

	}
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
		for (vector<segment>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			double currGam = segit->rateVec->back();
			if (maxGam > currGam)
			{
				segit->tbl->push_back(tableElem(1.0, 0.0, maxGam - currGam, i, j));
				segit->rateVec->push_back(maxGam);
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
double getRand()
{
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
@return The sum of all the rates calculated for the segment. For transition purposes
*/
double updateSegTable(shared_ptr<vector<CNT>> CNT_List, vector<segment>::iterator seg, double maxDist)
{
	double rate = 0;
	//iterate over CNTs
	int i = 0; //CNT index counter
	//originally structured without i and j
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		int j = 0; //segment index counter
		//iterate over all segments considered for seg
		for (vector<segment>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			double r;
			//Check if within range
			if ( ((r = tableElem::calcDist(seg->mid, segit->mid)) <= maxDist) && (r != 0))
			{
				auto theta = tableElem::calcThet(seg, segit);
				auto g = 6.4000e+19; //First draft estimate
				seg->tbl->push_back(tableElem(r,theta,g,i,j)); //tbl initialized in CNT::calculateSegments
				seg->rateVec->push_back(rate+=(seg->tbl->back()).getRate());//tbl initialized in CNT::calculateSegments
			}
			j++;
		}
		i++;
	}
	return rate;
}


/**
Converts numbers with some units to angstroms

@param unit The current unit
@param val The current value
@return the value in angstroms
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