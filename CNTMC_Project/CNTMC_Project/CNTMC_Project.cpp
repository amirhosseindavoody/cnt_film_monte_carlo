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


using namespace std;

//method declarations
string folderPathPrompt(bool incorrect);

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
	struct dirent *ent;
	list<string>*  fileList = new list<string>(0);
	int numFiles = 0;

	//Check if folder can be opened - should work due to above checks
	if ((resDir = opendir(resultFolderPath.c_str())) != nullptr)
	{
		//throw away first two results as they are . and ..
		readdir(resDir);readdir(resDir);
		//iterate over all of the real files
		while ((ent = readdir(resDir)) != nullptr)
		{
			fileList->push_back(ent->d_name);
			numFiles++;
		}
		closedir(resDir);
	} else
	{
		cout << "Could not open directory. Please try program again.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	//Iterate through the files and extract
	vector<CNT> *CNT_List = new vector<CNT>(numFiles);

	for (list<string>::iterator it = fileList->begin(); it != fileList->end(); ++it)
	{
		CNT_List->push_back(CNT(*(it), resultFolderPath));
	}
	

	return 0;
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
	delete inputResultsFolderPathArray;
	return returnString;
}