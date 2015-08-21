#pragma once
#include "tableElem.h"

struct energyRates
{
	energy src_energy;
	shared_ptr<vector<tableElem>> tbl;
	shared_ptr<vector<double>> rateVec;
};

