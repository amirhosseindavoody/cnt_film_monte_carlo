#ifndef file_management_h
#define file_management_h

#include <iostream>
#include <string>

using namespace std;

class file_management
{
	public:
		file_management();
		~file_management();
		
		string folderPathPrompt(bool incorrect);
		string checkPath(string path, bool folder);
		string xmlFilePathPrompt(bool incorrect);


		string resultFolderPath;
		string outputPath;
		
};


#endif // file_management_h