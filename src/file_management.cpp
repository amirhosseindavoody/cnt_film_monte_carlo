#include "stdafx.h"
#include <iostream>
#include <string>
#include <sys/stat.h>
#include "file_management.h"

using namespace std;

file_management:: file_management()
{
}

file_management:: ~file_management()
{
}

// Gets the path of the folder containing CNT mesh results
string file_management::folderPathPrompt(bool incorrect)
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
	char *inputResultsFolderPathArray = new char[inputPathLengthMax];
	cout << "Enter path of results folder:\n";
	if (incorrect)
		cin.ignore(); //if reentering, must ignore the next input
	cin.getline(inputResultsFolderPathArray, inputPathLengthMax);
	string returnString = inputResultsFolderPathArray;
	delete[] inputResultsFolderPathArray;
	return returnString;
}


// Checks to see if the file/folder exists and is accessable
string file_management::checkPath(string path, bool folder)
{

	//Checking validity of folder path
	while (true)
	{
		//Container for file/folder information
		struct stat info;

		if (stat(path.c_str(), &info) != 0)
		{
			printf("Cannot access %s\n", path.c_str());
			//get new path for next itr or quit
			if (folder){ path = folderPathPrompt(true); }
			else { path = xmlFilePathPrompt(true);}
		}
		//path is accessable
		else
		{
			//Check to see if correctly dir or file
			if ((folder && (info.st_mode & S_IFDIR)) || 
				(!folder && (info.st_mode & S_IFREG))) { return path; }
			
			//Incorrect type for the path. Reprompt or quit
			if (folder)
			{
				printf("%s is not a directory\n", path.c_str());
				path = folderPathPrompt(true);
			} 
			else
			{
				printf("%s is not a file\n", path.c_str());
				path = xmlFilePathPrompt(true);
			}
		}
	}
}


/**
Gets the path of the file used to configure the CNT mesh generation

@param incorrect Runs special prompt in case that cmd args were incorrect
@return The string containing the XML file path.
*/
string file_management::xmlFilePathPrompt(bool incorrect)
{
	if (incorrect)
	{
		string temp = " ";
		cout << "Re-enter XML config file path? [y/n]: ";
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
	cout << "Enter path of XML config file:\n";
	if (incorrect)
		cin.ignore(); //if reentering, must ignore the next input
	cin.getline(inputXMLFilePathArray, inputPathLengthMax);
	returnString = inputXMLFilePathArray;
	delete[] inputXMLFilePathArray;
	return returnString;
}