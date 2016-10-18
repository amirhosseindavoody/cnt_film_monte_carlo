#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>

using namespace std;

template <class T>
void write_log(T input)
{

	ofstream log_file;
	log_file.open("log.dat", ios::app);
	log_file << input << endl;

	if (log_file.fail())
	{
		cout << "error in writing to the log file!!!" << endl;
		exit(EXIT_FAILURE);
	}

	log_file.close();

	cout << input << endl;
}