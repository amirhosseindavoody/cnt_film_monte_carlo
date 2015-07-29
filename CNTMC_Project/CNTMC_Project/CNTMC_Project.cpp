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
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include <thread>
#include <omp.h>
#include <regex>


using namespace std;

//#define TFAC 2.302585092994046 //Factor to multiply 1/gamma by to get 90% of rand numbers giving tr less than deltaT
#define TFAC 1.203972804325936 //Factor to multiply 1/gamma by to get 70% of rand numbers giving tr less than deltaT

//method declarations
string folderPathPrompt(bool incorrect);
string xmlFilePathPrompt(bool incorrect);
string outputFolderPathPrompt(bool incorrect);
string checkPath(string path, bool folder);
double updateSegTable(shared_ptr<vector<CNT>> CNT_List, shared_ptr<segment> seg, 
	double maxDist, shared_ptr<vector<vector<int>>> heatMap, shared_ptr<vector<double>> rs, shared_ptr<vector<double>> thetas);
int getIndex(shared_ptr<vector<double>> vec, double val);
int getIndex(shared_ptr<vector<double>> vec, double val, int left, int right);
double getRand(bool excludeZero);
void addSelfScattering(shared_ptr<vector<CNT>> CNT_List, double maxGam);
void assignNextState(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> e, double gamma, shared_ptr<vector<double>> regionBdr);
double convertUnits(string unit, double val);
shared_ptr<vector<double>> linspace(double low, double high, int num);
void initRandomNumGen();
void injectExciton(shared_ptr<exciton> exciton, shared_ptr<vector<shared_ptr<segment>>> inContact);
bool hasMovedToOutContact(shared_ptr<exciton> exciton, shared_ptr<vector<double>> regionBdr, shared_ptr<vector<CNT>> CNT_List);
void markCurrentExcitonPosition(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> exciton, shared_ptr<vector<int>> currCount,
	shared_ptr<vector<double>> regionBdr);
void writeStateToFile(shared_ptr<ofstream> file, shared_ptr<vector<int>> currCount, double T);
void writeExcitonDistSupportingInfo(string outputPath, int numExcitons, double Tmax, double deltaT, double segLenMin, int numRegions,
	double xdim, double minBin, double rmax, int numBins, double lowAng, double highAng, int numAng, UINT64 numTSteps);
void updateExcitonList(int numExcitonsAtCont, shared_ptr<vector<shared_ptr<exciton>>> excitons, shared_ptr<vector<int>> currCount,
	shared_ptr<vector<shared_ptr<segment>>> inContact);

//Global variables
double ymax = 0; //stores maximum height the cylinders of the CNTs are found at. All will be greater than 0.

//Runs the file input, monte carlo, and file output sections of code
int main(int argc, char *argv[])
{
	unsigned int NUM_THREADS = omp_get_max_threads();

	//Varible initialization
	double segLenMin = 100.0; //[Angstroms]

	//Initialize random number generator before anything to ensure that getRand() always works
	initRandomNumGen();

	bool done = false; //Reused boolean variable for looping
	string resultFolderPath = " ";
	string inputXMLPath = " ";

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

	resultFolderPath = checkPath(resultFolderPath, true); //check result folder path


	////////////////////////////////// OUTPUT FOLDERS ///////////////////////////////

	//Outputfolder should be the same as the result folder for the mesh
	// All checks are already completed
	string outputPath = resultFolderPath + "/"; 

	//Results folder exists and can be accessed
	//grab file list
	DIR *resDir;
	struct dirent *ent = nullptr;
	unique_ptr<list<string>>  fileList(new list<string>(0));

	/*
	The next block creates a list of files in the directory chosen to be looked at.
	*/

	//////////////////////////// BUILD FILE LIST ///////////////////////////////////////////
	regex rgx("CNT_Num_\\d+\\.csv"); //files we want to look through
	//Check if folder can be opened - should work due to above checks
	if ((resDir = opendir(resultFolderPath.c_str())) != nullptr)
	{
		//throw away first two results as they are . and ..
		readdir(resDir); readdir(resDir);
		//iterate over all of the real files
		while ((ent = readdir(resDir)) != nullptr)
		{
			smatch matches; //match_results for string objects
			regex_search(string(ent->d_name), matches, rgx);
			if (!matches.empty())
			{
				fileList->push_back(ent->d_name);
			}
		}
		closedir(resDir); //deletes pointer
	}
	else
	{
		cout << "Could not open directory. Please try program again.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}
	delete ent;

	/////////////////////////////////// XML FILE PARSE /////////////////////////////////////
	
	//First need to get correct xml file name
	auto xmlNamePos = resultFolderPath.rfind('/', resultFolderPath.size() - 1) + 1;
	inputXMLPath = outputPath + resultFolderPath.substr(xmlNamePos, outputPath.size() - xmlNamePos) + ".xml";
	
	double rmax = 0; //maximum possible differences between sections of CNTs
	double xdim = 0; //Dimension in which exciton populations will be monitored
	while (!done)
	{
		try
		{
			rapidxml::xml_document<> doc; //create xml object
			rapidxml::file<> xmlFile(inputXMLPath.c_str()); //open file
			doc.parse<0>(xmlFile.data()); //parse contents of file
			rapidxml::xml_node<>* currNode = doc.first_node(); //gets the node "Document" or the root node
			currNode = currNode->first_node(); //Output folder
			//Speed up to device dimensions
			for (int i = 0; i < 6; i++)
			{
				currNode = currNode->next_sibling();
			}
			// DEVICE DIMENSIONS NODE //
			xdim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->value()));
			auto ydim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->value()));
			auto zdim = convertUnits(string(currNode->first_node()->value()),
				atof(currNode->first_node()->next_sibling()->next_sibling()->next_sibling()->value()));
			//incorrect units
			if (xdim == INT_MIN / 2.0 || ydim == INT_MIN / 2.0 || zdim == INT_MIN / 2.0)
			{
				printf("Configuration Error: Incorrect units for device dimensions.\nRefer to manual"
					" for valid unit entries.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			//incorrect range
			else if (xdim <= 0 || ydim <= 0 || zdim <= 0)
			{
				printf("Configuration Error: Must enter positive device dimensions.\n");
				system("pause");
				exit(EXIT_FAILURE);
			}
			// END DEVICE DIMENSIONS NODE //

			//x and z dim are the ground dimensions and y is height. Need to make sure bottom corner
			// to top other corner is included in r range. ymax is not defined so this is intermediate
			// value. It is changed after CNT initialization.
			rmax = sqrt(pow(xdim,2)+pow(ydim,2));

			// Will not let me delete, says it is null pointer even though I can get values from it. 
			//delete currNode; 
			done = true;
		}
		catch (runtime_error err)
		{
			string temp;
			int xmlArrayLength = 260; //Maximum path length for Windows
			cout << err.what();
			cout << "\n";
			cout << "Continue? [y/n]: ";
			cin >> temp;
			if (temp.compare("y") != 0)
			{
				system("pause");
				exit(EXIT_FAILURE);
			}
			char *inputXMLPathArray = new char[xmlArrayLength];
			system("cls");
			cout << "Enter config xml path (Example in program files directory):\n";
			cin.ignore();
			cin.getline(inputXMLPathArray, xmlArrayLength);
			inputXMLPath = inputXMLPathArray;
			delete inputXMLPathArray;
		}
	}


	/*
	This section creates the list of CNTs with all of the relevant information provided in their
	respective files getting processed.
	*/

	//////////////////////////// BUILD CNT AND SEGS //////////////////////////////////////////

	//Iterate through the files and extract
	shared_ptr<vector<CNT>> CNT_List(new vector<CNT>(0));

	for (list<string>::iterator it = fileList->begin(); it != fileList->end(); ++it)
	{
		CNT_List->push_back(CNT(*(it), resultFolderPath, segLenMin));
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

	//Set final rmax value
	rmax = sqrt(pow(rmax, 2) + pow(ymax,2));

	//////////////////////////// CREATE DIST VECTORS AND ARRAYS ///////////////////////////////

	shared_ptr<vector<double>> rs; //Vector containing range of r's
	rmax = 100 * ceil(rmax / 100.0);
	int numBins = static_cast<int>(rmax / 10.0); //number of bins to place r's into 
	double minBin = rmax / static_cast<double>(numBins); //[Angstroms] The size of the bins
	rs = linspace(minBin, rmax, numBins); //Builds rs vector within valid r range

	//builds angle vector from 1 to 90 degrees. Enough bins to cover all relevant angles
	// do not forget to use radians
	double lowAng = 1 * M_PI / 180.0;
	double highAng = 90 * M_PI / 180.0;
	int numAng = 90; //Number of angles to record. One per degree
	shared_ptr<vector<double>> thetas = linspace(lowAng, highAng, numAng);

	shared_ptr<vector<vector<int>>> heatMap(new vector<vector<int>>(rs->size()));
	//initialize all other vectors in the vector to the correct size
	for (vector<vector<int>>::iterator it = heatMap->begin(); it != heatMap->end(); ++it)
	{
		it->resize(thetas->size());
	}

	/*
	Building the table: Each CNT has a list of segments which has an empty list of table element objects.
	The next section will be filling the segments vector of table elements by calculating the necessary
	additions to it.
	*/

	//////////////////////////// BUILD TABLE ////////////////////////////////////////////////

	///////////// Variables for placing excitons in the future /////////////
	int numRegions = static_cast<int>(xdim / segLenMin); //number of regions in the simulation
	double extra = xdim - numRegions*segLenMin; //extra amount * numRegions that should be added to segLenMin
	double regLen = segLenMin + extra / numRegions; //The length of each region.

	//The boundary of the rgions in the x direction
	shared_ptr<vector<double>> regionBdr = linspace(regLen - (xdim / 2), xdim / 2, numRegions);
	shared_ptr<vector<int>> secCountPerReg = make_shared<vector<int>>(vector<int>(numRegions)); //To get dist stats
	//List of segments in the first region, which is used as a input contact
	shared_ptr<vector<shared_ptr<segment>>> inContact(new vector<shared_ptr<segment>>(0));
	//List of segments in the last region, which is used as an exit contact
	shared_ptr<vector<shared_ptr<segment>>> outContact(new vector<shared_ptr<segment>>(0));

	
	//iterate through all of the CNTs and segments
	double maxDist = 300; //[Angstroms]
	//The total number of segments in the simulation, used in exciton placement
	auto numSegs = 0;
	//The maximum of sums of gammas from each segment. This sets the constant gamma value for the entire simulation
	double gamma = 0;

	omp_set_num_threads(NUM_THREADS);
	int nThreads = 0;

	//loop over CNTs
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		//define variables before using in parallel processing loop
		int begin = 0;
		int end = static_cast<int>(cntit->segs->size());
		int segidx;
		shared_ptr<segment> currSeg;
		double newGamma;
		int regIdx;

		#pragma omp parallel default(shared) private(segidx,newGamma,regIdx,currSeg) shared(gamma, begin, end,numSegs)
		{
			#pragma omp master
			nThreads = omp_get_num_threads();
			cout << nThreads << endl;
			
			//loop over segments in each CNTs
			#pragma omp parallel for
			for (segidx = begin; segidx < end; segidx++)
			{
				currSeg = (*(cntit->segs))[segidx];
				regIdx = getIndex(regionBdr, currSeg->mid(0));
				(*secCountPerReg)[regIdx]++; //increment the count based on where section is
				if (regIdx == 0){ inContact->push_back(currSeg); } //First region is injection contact
				else if (regIdx == secCountPerReg->size() - 1){ outContact->push_back(currSeg); } //last region is output contact
				//get add to each segment relevant table entries
				newGamma = updateSegTable(CNT_List, currSeg, maxDist, heatMap, rs, thetas);
				#pragma omp critical
				if (newGamma > gamma)
				{
					gamma = newGamma;
				}
				#pragma omp atomic
				numSegs++;
			}
		}
		if (nThreads == NUM_THREADS) {
			cout << "Tube Number: " << cntit->getCNTNum() << endl;
			printf_s("%d OpenMP threads were used.\n", NUM_THREADS);
		}
	}

	/////////////////////////////// OUTPUT HEATMAP //////////////////////////////////////

	string heatMapFileName = outputPath + "heatMap.csv";
	ofstream heatMapFile;
	heatMapFile.open(heatMapFileName);
	//iterate through all of the rs and then thetas while printing to file
	for (UINT32 i = 0; i < rs->size(); i++)
	{
		for (UINT32 j = 0; j < thetas->size(); j++)
		{
			heatMapFile << (*heatMap)[i][j] << ",";
		}
		heatMapFile << "\n";
	}
	heatMapFile.close();



	/////////////////////////////// ADD SELF SCATTERING //////////////////////////////////

	addSelfScattering(CNT_List, gamma);

	/*
	Now that the tables have been built, the next step is to populate the mesh with excitons. The way this will happen
	is, after a specific number of excitons are chosen to be created, each exciton will be created and assigned an index
	that corresponds to the CNT mesh. The available locations for excitons to start will be determined by the number of 
	segments that are within the delta x region specified. The exciton assignment will be equally and randomly distributed 
	among the list of segments in the delta x region.
	*/

	//////////////////////////// PLACE EXCITONS ////////////////////////////////////////////

	//The number of excitons to be in injection contact
	int numExcitonsAtCont = 10000;
	//Vector of excitons. Positions and energies must still be assigned
	shared_ptr<vector<shared_ptr<exciton>>> excitons(new vector<shared_ptr<exciton>>(numExcitonsAtCont));
	for (int exNum = 0; exNum < numExcitonsAtCont; exNum++)
	{
		(*excitons)[exNum] = make_shared<exciton>(exciton()); //initialize exciton at location in exciton list
		injectExciton((*excitons)[exNum], inContact);
	}


	//////////////////////////////////// TIME STEPS ///////////////////////////////////
	double deltaT = (1 / gamma)*TFAC; //time steps at which statistics are calculated
	double Tmax = deltaT * 100000; //maximum simulation time
	double T = 0; //Current simulation time, also time at which next stats will be calculated

	shared_ptr<vector<int>> currCount; //The count of excitons in each region of the CNT mesh

	//File output initializations
	string excitonDistFileName = outputPath + "excitonDist.csv";
	shared_ptr<ofstream> excitonDistFile(new ofstream);
	excitonDistFile->open(excitonDistFileName);

	
	/*
	This section will consist of iterating until the maximum time has been reached. Each iteration
	for T will contain a single/multiple step for each exciton. 
	*/
	while (T <= Tmax) //iterates over time
	{
		//reset the exciton count after each time step
		currCount = make_shared<vector<int>>(vector<int>(numRegions));
		T += deltaT; //set new time checkpoint
		for (UINT32 exNum = 0; exNum < excitons->size(); exNum++) //iterates over excitons once
		{
			/*There is a change that the previous tr that was calculated was so long that it
			not only has extra time in the next deltaT, but it skips it completely. In this case
			the exciton movement is skipped all together and its extra time is decreased by deltaT
			until the exciton can move again. Otherwise the code runs as usual.
			*/
			double extraT = (*excitons)[exNum]->getTExtra();
			if (extraT > deltaT)
			{
				(*excitons)[exNum]->setTExtra((*excitons)[exNum]->getTExtra() - deltaT);
				markCurrentExcitonPosition(CNT_List, (*excitons)[exNum], currCount, regionBdr);
			}
			else
			{
				double tr_tot = extraT; //the sum of all tr's in the current deltaT time step
				/*give time to excitons that have no extra time. This happens only when excitons
				are just injected.*/
				if (extraT == 0)
				{
					tr_tot += -(1 / gamma)*log(getRand(true));
				}
				while (tr_tot <= deltaT) 
				{
					//choose new state
					assignNextState(CNT_List, (*excitons)[exNum], gamma, regionBdr);
					tr_tot += -(1 / gamma)*log(getRand(true)); // add the tr calculation to current time for individual particle
				}
				//Recording distribution
				markCurrentExcitonPosition(CNT_List, (*excitons)[exNum], currCount, regionBdr);
				//record the time past deltaT that the assign next state will cover
				(*excitons)[exNum]->setTExtra(tr_tot - deltaT);
			}
		}
		//Output count vector to file since we want results after each time step.
		writeStateToFile(excitonDistFile, currCount, T);

		//Update Exciton List for injection and exit contact
		updateExcitonList(numExcitonsAtCont, excitons, currCount, inContact);
	}
	//Close files finish program
	excitonDistFile->close();

	UINT64 numTSteps = static_cast<UINT64>(T / deltaT);

	//Write helper information
	writeExcitonDistSupportingInfo(outputPath, numExcitonsAtCont, Tmax, deltaT, segLenMin, numRegions, xdim,
		minBin, rmax, numBins, lowAng, highAng, numAng,numTSteps);

	return 0;
}

/**
Removes excitons from the list if they are in the exit contact and inject excitons into the inContact if there are
not enough excitons

@param numExcitonAtCont The number of excitons at the injection contact
@param excitons The list of excitons to modify
@param currCount The count of excitons in each of the regions
@param inContact The list of segments in the injection contact
*/
void updateExcitonList(int numExcitonsAtCont, shared_ptr<vector<shared_ptr<exciton>>> excitons, 
	shared_ptr<vector<int>> currCount, shared_ptr<vector<shared_ptr<segment>>> inContact)
{
	// Adding excitons to the injection contact
	int numExAdd = 0; //The number of excitons to be added to the injection contact
	if ((numExAdd = numExcitonsAtCont - (*currCount)[0]) > 0)
	{
		for (int i = 0; i < numExAdd; i++)
		{
			excitons->push_back(make_shared<exciton>(exciton())); //initialize exciton at end of exciton list
			injectExciton((*excitons)[excitons->size()-1], inContact); //last element is the exciton to inject
		}
	}

	//Removing excitons from the exit contact
	int numExRem = (*currCount)[currCount->size() - 1]; //Last element in the count list = num of excitons to remove
	if (numExRem > 0)
	{
		for (int i = 0; i < excitons->size(); i++)
		{
			if ((*excitons)[i]->isAtOutContact())
			{
				shared_ptr<exciton> swap = (*excitons)[excitons->size() - 1];
				(*excitons)[excitons->size() - 1] = (*excitons)[i];
				(*excitons)[i] = swap;
				excitons->pop_back();
				//excitons->erase(excitons->begin() + i);
				i--;
			}
		}
	}
}





/**
Writes all passed information to a file that can be read by matlab to get the appropriate variable declarations
*/
void writeExcitonDistSupportingInfo(string outputPath, int numExcitons, double Tmax, double deltaT, double segLenMin, int numRegions, 
	double xdim, double minBin, double rmax, int numBins, double lowAng, double highAng, int numAng, UINT64 numTSteps)
{
	string detailsFileName = outputPath + "details.csv";
	shared_ptr<ofstream> detailsFile(new ofstream);
	detailsFile->open(detailsFileName);
	*detailsFile << "numExcitons,Tmax,deltaT,segLenMin,numRegions,xdim,minBin,rmax,numBins,lowAng,highAng,numAng,numTSteps" << endl;
	*detailsFile << numExcitons << "," << Tmax << "," << deltaT << "," << segLenMin << "," << numRegions << "," << xdim 
		<< "," << minBin << "," << rmax << "," << numBins << "," << lowAng << "," << highAng << "," << numAng << 
		"," << numTSteps << endl;
}

/**
Writes the pertinent information of the current state of simulation to file

@param file The file to write to
@param currCount The distribution information to write to file
@param T The current time of the simulation
*/
void writeStateToFile(shared_ptr<ofstream> file, shared_ptr<vector<int>> currCount, double T)
{
	*file << T << "," << (*currCount)[0]; //Time followed by the counts
	for (UINT32 i = 1; i < currCount->size(); i++)
	{
		*file << "," << (*currCount)[i];
	}
	*file << endl;
}


/**
Takes the current exciton and adds it to the count of excitons in a certain region. 
This is called once per exciton per time step.

@param CNT_List The indexable list that contains the position information of exciton
@param exciton The exciton that we want to mark the position of.
@param currCount The vector that stores the count of excitons in the regions
@param regionBdr The region vector that allows selection of currCount index
*/
void markCurrentExcitonPosition(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> exciton, 
	shared_ptr<vector<int>> currCount, shared_ptr<vector<double>> regionBdr)
{
	//add count to the currCount vector in the location corresponding to the exciton 
	//   region location
	(*currCount)[getIndex(regionBdr,
		((*(*CNT_List)[exciton->getCNTidx()].segs)[exciton->getSegidx()]->mid(0)))]++;
}

/**
Checks to see if an exciton has moved into the out contact

@param exciton The exciton to check
@param regionBdr The array that stores contact information
@param CNT_List The indexable list that will reveal exciton location
@return True if the exciton has moved into the output contact, false otherwise
*/
bool hasMovedToOutContact(shared_ptr<exciton> exciton, shared_ptr<vector<double>> regionBdr, 
	shared_ptr<vector<CNT>> CNT_List)
{
	//If x component is creater than the second to last index of regionBdr, then it is in the output
	return (((*(*CNT_List)[exciton->getCNTidx()].segs)[exciton->getSegidx()]->mid(0)) 
					>= ((*regionBdr)[regionBdr->size() - 2]));
}


/**
Places the specified exciton into the input contact

@param exciton The exciton to add to the contact
@param inContact The list of segments for the injection contact
*/
void injectExciton(shared_ptr<exciton> exciton, shared_ptr<vector<shared_ptr<segment>>> inContact)
{
	exciton->setEnergy(static_cast<int>(round(getRand(false)) + 1)); //randomly set the energy of the exciton
	//choose a destination segment
	shared_ptr<segment> injectedSeg = (*inContact)[static_cast<int>(rand() % inContact->size())];
	//The self scattering table will have the correct indices for the current segment
	shared_ptr<vector<tableElem>> tbl = injectedSeg->tbl;
	//Set the exciton indices to the current segment
	exciton->setCNTidx((*tbl)[tbl->size() - 1].getTubeidx());
	exciton->setSegidx((*tbl)[tbl->size() - 1].getSegidx());
	exciton->setAtOutContact(false); //initialized out contact boolean

}


/**
Assigns the specified exciton to the next state in the simulation.

@param CNT_List The list of carbon nanotubes
@param e The exciton to be updated
@param inContact The list of segments for the injection contact
@param regionBdr The array that will help decide whether or not the exciton is in the out contact
*/
void assignNextState(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> e, double gamma, shared_ptr<vector<double>> regionBdr)
{
	//Segment the current exciton is located on
	shared_ptr<segment> seg = (*((*CNT_List)[e->getCNTidx()].segs))[e->getSegidx()];
	//Get the table index for the
	int tblIdx = getIndex(seg->rateVec, getRand(false)*gamma);
	//stores information about the excitons destination
	tableElem tbl = (*seg->tbl)[tblIdx];
	/*
	It was decided that there are no limits on the number of excitons that
	can be on a segment. 7/20/15
	*/
	e->setCNTidx(tbl.getTubeidx());
	e->setSegidx(tbl.getSegidx());
	e->setAtOutContact(hasMovedToOutContact(e, regionBdr, CNT_List));

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
		for (vector<shared_ptr<segment>>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			double currGam = (*segit)->rateVec->back();
			if (maxGam > currGam)
			{
				(*segit)->tbl->push_back(tableElem(1.0, 0.0, maxGam - currGam, i, j));
				(*segit)->rateVec->push_back(maxGam);
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
double getRand(bool excludeZero)
{
	if (excludeZero)
	{
		int r;
		while ((r = rand()) == 0){}
		return static_cast<double>(r) / static_cast<double>(RAND_MAX);
	}
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
	return getIndex(vec, val, 0, static_cast<int>(vec->size() - 1));
	//size_t is 64 bit unsigned. int is 32 bit signed. Cast should still work always
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
@param colorMap A count of all rs at particular thetas to get mesh statistics
@param rs The r values that are needed to place values in the heat map
@param thetas The angles that are needed to place values in the heat map
@return The sum of all the rates calculated for the segment. For transition purposes
*/
double updateSegTable(shared_ptr<vector<CNT>> CNT_List, shared_ptr<segment> seg,
	double maxDist, shared_ptr<vector<vector<int>>> heatMap, shared_ptr<vector<double>> rs, shared_ptr<vector<double>> thetas)
{
	double rate = 0;
	//iterate over CNTs
	int i = 0; //CNT index counter
	//originally structured without i and j
	for (vector<CNT>::iterator cntit = CNT_List->begin(); cntit != CNT_List->end(); ++cntit)
	{
		int j = 0; //segment index counter
		//iterate over all segments considered for seg
		for (vector<shared_ptr<segment>>::iterator segit = cntit->segs->begin(); segit != cntit->segs->end(); ++segit)
		{
			double r = tableElem::calcDist(seg->mid, (*segit)->mid);
			auto theta = tableElem::calcThet(seg, segit);

			//Heat Map Additions
			if (r != 0)
			{
				(*heatMap)[getIndex(rs, r)][getIndex(thetas, theta)]++; //////// Heat Map ///////
				//Check if within range
				if (r <= maxDist) /////// Building TABLE /////
				{
					auto g = 6.4000e+19; //First draft estimate
					seg->tbl->push_back(tableElem(r, theta, g, i, j)); //tbl initialized in CNT::calculateSegments
					seg->rateVec->push_back(rate += (seg->tbl->back()).getRate());//tbl initialized in CNT::calculateSegments
				}
			}	
			j++;
		}
		i++;
	}
	return rate;
}


/**
Gets the path of the folder containing CNT mesh results

@param incorrect Runs special prompt in case that cmd args were incorrect
@return The string containing the results folder path.
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

/**
Gets the path of the file used to configure the CNT mesh generation

@param incorrect Runs special prompt in case that cmd args were incorrect
@return The string containing the XML file path.
*/
string xmlFilePathPrompt(bool incorrect)
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

/**
Gets the path of output folder

@param incorrect Runs special prompt in case that cmd args were incorrect
@return The string containing the XML file path.
*/
string outputFolderPathPrompt(bool incorrect)
{
	if (incorrect)
	{
		string temp = " ";
		cout << "Re-enter output folder path? [y/n]: ";
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
	cout << "Enter path of output folder:\n";
	if (incorrect)
		cin.ignore(); //if reentering, must ignore the next input
	cin.getline(inputXMLFilePathArray, inputPathLengthMax);
	returnString = inputXMLFilePathArray;
	delete[] inputXMLFilePathArray;
	return returnString;
}

/**
Converts numbers with some units to angstroms

@param unit The current unit
@param val The current value
@return the value in angstroms
*/
double convertUnits(string unit, double val)
{
	if (unit.compare("mm") == 0 || unit.compare("millimeter") == 0)
	{
		return val * 10000000;
	}
	else if (unit.compare("um") == 0 || unit.compare("micrometer") == 0)
	{
		return val * 10000;
	}
	else if (unit.compare("nm") == 0 || unit.compare("nanometer") == 0)
	{
		return val * 10;
	}
	else if (unit.compare("pm") == 0 || unit.compare("picometer") == 0)
	{
		return val * .01;
	}
	else if (unit.compare("A") == 0 || unit.compare("angstrom") == 0)
	{
		return val;
	}
	else
	{
		return INT_MIN;
	}
}

/**
Checks to see if the file/folder exists and is accessable

@param path The path of the file/folder
@param folder If true, then the path is for a folder, file otherwise
@return The string of the path that works
*/
string checkPath(string path, bool folder)
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
Creates a vector with matlab style linspace numbering

@param low The lowest number to be in vector
@param high The highest number to be in vector
@param num The number of points to be in vector
@return A pointer to the resulting vector
*/
shared_ptr<vector<double>> linspace(double low, double high, int num)
{
	shared_ptr<vector<double>> retVec(new vector<double>(num));
	double step = (high - low) / static_cast<double>(num - 1);
	for (int i = 0; i < num; i++)
	{
		(*retVec)[i] = low;
		low += step;
	}
	return retVec;
}

/**
Initialized the random number generator by providing a seed value from 
the time.
*/
void initRandomNumGen()
{
	//Initialize random number generation
	time_t seconds;
	time(&seconds); //assign time from clock
	//Seed the random number generator
	srand(static_cast<int>(seconds));
}