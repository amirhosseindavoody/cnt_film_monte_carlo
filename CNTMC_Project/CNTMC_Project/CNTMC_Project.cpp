// CNTMC_Project.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <string>
#include <sys/stat.h>


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