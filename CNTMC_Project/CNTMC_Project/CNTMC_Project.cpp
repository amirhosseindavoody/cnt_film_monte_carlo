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


using namespace std;

//method declarations
string folderPathPrompt(bool incorrect);
void updateSegTable(shared_ptr<vector<CNT>> CNT_List, vector<segment>::iterator seg, double maxDist);

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
	vector<double> cntProb = vector<double>(CNT_List->size());

	//Must set first index so that rest of the table can be built with simple loop
	vector<CNT>::iterator it = CNT_List->begin();
	cntProb[0] = it->segs->size() / static_cast<double>(numSegs);
	++it;

	auto cntIdx = 1; //keeps index of cntProb vector as iterate through CNT_List
	for ( it; it != CNT_List->end(); ++it)
	{
		cntProb[cntIdx] = it->segs->size() / static_cast<double>(numSegs) + cntProb[cntIdx-1];
		cntIdx++;
	}




	return 0;
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