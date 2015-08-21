#pragma once
#include "tableElem.h"
#include <vector>
#include <memory>

class tableElem;

using namespace std;

struct energyRates
{
	energy src_energy;
	shared_ptr<vector<tableElem>> tbl;
	shared_ptr<vector<double>> rateVec;
};
