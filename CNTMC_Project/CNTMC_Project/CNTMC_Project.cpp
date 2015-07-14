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

//method declarations
string folderPathPrompt(bool incorrect);
void updateSegTable(shared_ptr<vector<CNT>> CNT_List, vector<segment>::iterator seg, double maxDist);
int getIndex(shared_ptr<vector<double>> vec, double prob);
int getIndex(shared_ptr<vector<double>> vec, double prob, int left, int right);
double getRand();

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
	
	//iterate through all of the CNTs and segments
	double maxDist = 500; //[Angstroms]
	//The total number of segments in the simulation, used in exciton placement
	auto numSegs = 0;
	//loop over CNTs
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		//loop over segments in each CNTs
		for (vector<segment>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			//get add to each segment relelant table entries
			updateSegTable(CNT_List, segit, maxDist);
			numSegs++;
		}
	}

	/*
	Now that the tables have been built, the next step is to populate the mesh with excitons. The way this will happen
	is, after a specific number of excitons are chosen to be created, each exciton will be created and assigned an index
	that corresponds to the CNT mesh. First a random CNT will be chosen. This choice will come from a weighted probability
	distribution, the weight being placed on the number of segments each CNT has vs the total number of segments. After the
	CNT number is chosen, then the segment will be chosen from the segments part of the CNT.
	*/

	//Initialize random number generation
	time_t seconds;
	time(&seconds); //assign time from clock
	//Seed the random number generator
	srand(static_cast<int>(seconds));

	//Build the probability vector for CNTs
	shared_ptr<vector<double>> cntProb(new vector<double>(CNT_List->size()));

	//Must set first index so that rest of the table can be built with simple loop
	vector<CNT>::iterator it = CNT_List->begin();
	cntProb[0] = it->segs->size() / static_cast<double>(numSegs);
	++it;

	auto cntIdx = 1; //keeps index of cntProb vector as iterate through CNT_List
	for ( it; it != CNT_List->end(); ++it)
	{
		(*cntProb)[cntIdx] = it->segs->size() / static_cast<double>(numSegs) + cntProb[cntIdx-1];
		cntIdx++;
	}

	int numExcitons = 10; //number of excitons to add to simulation

	if (numExcitons > 2*numSegs)
	{
		cout << "Error: More excitons than possible number of places for them in simulation.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	shared_ptr<vector<exciton>> excitons(new vector<exciton>(numExcitons)); //vector that stores all excitons

	for (exciton &e : *excitons)
	{
		done = false; // Flag for successful exciton assignment
		//randomly sets energy level to 1 or 2
		e.setEnergy(round(getRand())+1); 

		//go until suitable position for exciton is found
		while (!done)
		{
			e.setCNTidx(getIndex(cntProb, getRand()));
		}
	}

	return 0;
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
Finds the index of the vector that has the number closest to but greater than prob.
Does this recursively.

@param vec The vector to search
@param prob The number to compare the indecies to
@return the index that requires the conditions in the method description
*/
int getIndex(shared_ptr<vector<double>> vec, double prob)
{	
	//parameter check
	if (vec == nullptr)
	{
		cout << "Error: Empty vector passed to getIndex()";
		system("pause");
		exit(EXIT_FAILURE);
	} else if (prob < 0 || prob > 1)
	{
		cout << "Error: Invalid probability range passed to getIndex().\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	//small simple cases that do not work with recursive helper method
	if (vec->size() == 1) {return 0;}
	if (vec->size() == 2)
	{
		if ((*vec)[0] >= prob){ return 0; }
		return 1;
	}

	//use recursive helper method
	return getIndex(vec, prob, 0, vec->size() - 1);
}

/**
The helper method of the 2 parameter get index method. Simple binary search

@param vec The vector to search
@param prob The number to compare the indecies to
@param right The highest index that is part of the vector
@param lef The lowest index that is part of the vector
@return the index that requires the conditions in the method description. Returns -1 if it fails
*/
int getIndex(shared_ptr<vector<double>> vec, double prob, int left, int right)
{
	if (right < left) { return left; }

	if (right == left)
	{
		if ((*vec)[right] < prob) { return right + 1; }
		return right; 
	} //base case
	
	int center = static_cast<int>(floor(static_cast<double>(right + left) / 2.0));
	if ((*vec)[center] == prob){ return center; } //base case
	
	//recursive cases
	if ((*vec)[center] > prob){ return getIndex(vec, prob, left, center - 1); }
	if ((*vec)[center] < prob){ return getIndex(vec, prob, center + 1, right); }
	return -1;

}

void updateSegTable(shared_ptr<vector<CNT>> CNT_List, vector<segment>::iterator seg, double maxDist)
{
	//iterate over CNTs
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		//iterate over all segments considered for seg
		for (vector<segment>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			double r;
			//Check if within range
			if ( (r = tableElem::calcDist(seg->mid, segit->mid)) <= maxDist)
			{
				auto theta = tableElem::calcThet(seg, segit);
				auto g = 6.4000e+19; //First draft estimate
				seg->tbl->push_back(tableElem(r,theta,g,cntit->getCNTNum()-1,segit->segNum));
			}
		}
	}
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