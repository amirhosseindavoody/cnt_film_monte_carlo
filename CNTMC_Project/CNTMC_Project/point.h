/**
point.h
Purpose: Struct to store data about a point

@author Alex Gabourie
@version 1.00
*/

#pragma once
#include <memory>

struct point
{
	//information for row
	uint32_t r_idx;
	double r_val;

	//information for column
	uint32_t t_idx;
	double t_val;

	//value gathered from index
	double val;

};