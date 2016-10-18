#ifndef functions_h
#define functions_h

#include <stdio.h>
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
#include <regex>
#include <omp.h>
#include <stdint.h>


//method declarations
double updateSegTable(shared_ptr<vector<CNT>> CNT_List, vector<shared_ptr<segment>>::iterator seg, 
	double maxDist, shared_ptr<vector<vector<int>>> heatMap, shared_ptr<vector<double>> rs, shared_ptr<vector<double>> thetas);
int getIndex(shared_ptr<vector<double>> vec, double val);
double getRand(bool excludeZero);
void addSelfScattering(shared_ptr<vector<CNT>> CNT_List, double maxGam);
void assignNextState(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> e, double gamma, shared_ptr<vector<double>> regionBdr);
double convert_units(string unit, double val);
shared_ptr<vector<double>> linspace(double low, double high, int num);
void init_random_number_generator();
void injectExciton(shared_ptr<exciton> exciton, shared_ptr<vector<shared_ptr<segment>>> inContact);
bool hasMovedToOutContact(shared_ptr<exciton> exciton, shared_ptr<vector<double>> regionBdr, shared_ptr<vector<CNT>> CNT_List);
void markCurrentExcitonPosition(shared_ptr<vector<CNT>> CNT_List, shared_ptr<exciton> exciton, shared_ptr<vector<int>> currCount,
	shared_ptr<vector<double>> regionBdr);
void writeStateToFile(shared_ptr<ofstream> file, shared_ptr<vector<int>> currCount, double T);
void writeExcitonDistSupportingInfo(string outputPath, int numExcitons, double Tmax, double deltaT, double segLenMin, int numRegions,
	double xdim, double minBin, double rmax, int numBins, double lowAng, double highAng, int numAng, uint64_t numTSteps, double regLenMin, string runtime);
void updateExcitonList(int numExcitonsAtCont, shared_ptr<vector<shared_ptr<exciton>>> excitons, shared_ptr<vector<int>> currCount,
	shared_ptr<vector<shared_ptr<segment>>> inContact);


#endif // functions_h