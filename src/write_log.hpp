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

	const int max_path_length = 256;

	char buffer[max_path_length];

	char *path = getcwd(buffer, max_path_length);
	string current_path = path;
	cout << current_path << endl;

}