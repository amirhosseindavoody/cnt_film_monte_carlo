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


using namespace std;

//method declarations
string folderPathPrompt(bool incorrect);
void updateSegTable(shared_ptr<vector<CNT>> CNT_List, shared_ptr<vector<tableElem>> tblEl);

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
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		for (vector<segment>::iterator segit = cntit->segs.begin(); segit != cntit->segs.end(); ++segit)
		{
			updateSegTable(CNT_List, segit->tbl);
		}
	}

	return 0;
}

void updateSegTable(shared_ptr<vector<CNT>> CNT_List, shared_ptr<vector<tableElem>> tblEl)
{
	string yo = "dude man";
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