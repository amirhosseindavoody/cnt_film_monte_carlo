#include <iostream>
#include <limits>

#include "ff.h"

namespace mc
{

// // constructor
// free_flight::free_flight(mc::arr1d accel)
// {
// 	// make sure that there is an inverse for acceleration value. othersize put a small value instead of acceleration so that 1/accel is not infinity
// 	for (auto &elem : accel)
// 	{
// 		if (! std::isfinite(1./elem))
// 		{
// 			elem = std::numeric_limits<mc::t_float>::epsilon();
// 		}
// 	}
//
// 	for (int i=0; i<accel.size(); ++i)
// 	{
// 		_acceleration[i] = accel[i];
// 		_inv_acceleration[i] = 1./accel[i];
// 	}
// };



} // mc namespace
