#include <iostream>
#include <memory>
#include <experimental/filesystem>
#include <ctime>
#include <cmath>
#include <array>

#include "utility.h"

namespace mc
{

// initialize the seed of the random number generator
void init_random_number_generator()
{
	std::srand(std::time(0));
};

} // namespace mc