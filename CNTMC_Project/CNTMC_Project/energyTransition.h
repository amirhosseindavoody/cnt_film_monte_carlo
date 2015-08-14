#pragma once
#include <memory>
#include <vector>

using namespace std;

//Type to store energy information
enum energy
{
	E11,
	E22
};

class energyTransition
{

	energy source_energy; //energy of exciton on source segment
	energy dest_energy; //energy of exciton on destination segment
	shared_ptr<vector<shared_ptr<vector<double>>>> rtTable; // r and theta table values

public:
	energyTransition();
	~energyTransition();
};

