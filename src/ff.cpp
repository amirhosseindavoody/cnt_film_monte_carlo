#include <iostream>
#include <limits>

#include "ff.h"

namespace mc
{

// constructor
free_flight::free_flight(mc::arr1d accel)
{
	// make sure that there is an inverse for acceleration value. othersize put a small value instead of acceleration so that 1/accel is not infinity
	for (auto &elem : accel)
	{
		if (! std::isfinite(1./elem))
		{
			elem = std::numeric_limits<mc::t_float>::epsilon();
			// elem = 1.0;
		}
	}

	for (int i=0; i<accel.size(); ++i)
	{
		_acceleration[i] = accel[i];
		_inv_acceleration[i] = 1./accel[i];
	}

	std::cout << "\n****\n";
	std::cout << "acceleration = ";
	for (auto elem : _acceleration)
		std::cout << elem << " , ";
	std::cout << std::endl;

	std::cout << "inverse acceleration = ";
	for (auto elem : _inv_acceleration)
		std::cout << elem << " , ";
	std::cout << std::endl;
	std::cout << "****\n";

	_number_of_flights_without_boundary_collisions = 0;
};

// perform free flight
void free_flight::fly(mc::arr1d &pos, mc::arr1d &velocity, const mc::t_float &eff_mass, const mc::t_float &dt, const mc::arr1d& volume)
{
	const mc::t_float coeff = dt*dt/2.;

	for (int i=0; i<pos.size(); ++i)
	{
		pos[i] +=  velocity[i]*dt + _acceleration[i]*coeff;
		velocity[i] += _acceleration[i]*dt;
	}
};

// get constant reference to acceleration vector
const mc::arr1d& free_flight::acceleration()
{
	return _acceleration;
};

// check for collision to boundaries
void free_flight::check_boundary(mc::arr1d &pos, mc::arr1d &velocity, const mc::arr1d& old_pos, const mc::arr1d& old_velocity, const mc::t_float &eff_mass, const mc::t_float &dt, const mc::arr1d& volume)
{
	mc::t_float t_collision;
	mc::t_float temp_expr;

	bool collision;
	mc::t_float o_pos, t_pos, n_pos;
	mc::t_float o_vel, t_vel, n_vel;
	mc::t_float expression2, expression1;

	_number_of_flights_without_boundary_collisions += 1;


	// for (int i=0; i<pos.size(); ++i)
	for (int i=0; i<1; ++i)
	{
		collision = false;
		o_pos = old_pos[i];
		t_pos = pos[i];
		o_vel = old_velocity[i];
		t_vel = velocity[i];

		if (pos[i] > volume[i])
		{
			temp_expr = 2*_acceleration[i]*(old_pos[i]-volume[i])/std::pow(old_velocity[i],2);
			if (temp_expr < 1.e-3)
			{
				
			}
			t_collision = (-old_velocity[i]+std::sqrt(std::pow(old_velocity[i],2)-2.*_acceleration[i]*(old_pos[i]-volume[i])))*(_inv_acceleration[i]);
			pos[i] = volume[i] + (-old_velocity[i]-_acceleration[i]*t_collision)*(dt-t_collision) + _acceleration[i]*std::pow(dt-t_collision,2)/2.;
			velocity[i] = ((-old_velocity[i]-_acceleration[i]*t_collision)+_acceleration[i]*(dt-t_collision));

			// under_sqrt = std::sqrt(std::pow(old_velocity[i],2)-2.*_acceleration[i]*(old_pos[i]-volume[i]));
			expression1 = std::pow(old_velocity[i],2);
			expression2 = 2.*_acceleration[i]*(old_pos[i]-volume[i]);
			collision = true;
		}
		else if(pos[i] < 0)
		{
			t_collision = (-old_velocity[i]-std::sqrt(std::pow(old_velocity[i],2)-2.*_acceleration[i]*(old_pos[i])))*(_inv_acceleration[i]);
			pos[i] = (-old_velocity[i]-_acceleration[i]*t_collision)*(dt-t_collision) + _acceleration[i]*std::pow(dt-t_collision,2)/2.;
			velocity[i] = ((-old_velocity[i]-_acceleration[i]*t_collision)+_acceleration[i]*(dt-t_collision));

			// under_sqrt = std::sqrt(std::pow(old_velocity[i],2)-2.*_acceleration[i]*(old_pos[i]));
			expression1 = std::pow(old_velocity[i],2);
			expression2 = 2.*_acceleration[i]*(old_pos[i]);

			collision = true;
		}

		n_pos = pos[i];
		n_vel = velocity[i];

		if (collision)
		{
			std::cout << std::scientific;
			std::cout << "***\nboundary collision detected:\n";
			std::cout << "number of steps without boundary collision= " << _number_of_flights_without_boundary_collisions << '\n';
			std::cout << "dim= " << i << std::endl;
			std::cout << " ,old position= " << o_pos << std::endl;
			std::cout << " ,tem position= " << t_pos << std::endl;
			std::cout << " ,new position= " << n_pos << std::endl;
			std::cout << "\n";
			std::cout << " ,old velocity= " << o_vel << std::endl;
			std::cout << " ,tem velocity= " << t_vel << std::endl;
			std::cout << " ,new velocity= " << n_vel << std::endl;
			std::cout << "\n";
			std::cout << " ,expression1= " << expression1 << std::endl;
			std::cout << " ,expression2= " << expression2 << std::endl;
			std::cout << " ,t_collision= " << t_collision << "\n";
			std::cout << "\nlimits:" << '\n';
			std::cout << "   ,min: " << std::numeric_limits<mc::t_float>::min() << '\n';
			std::cout << "   ,max: " << std::numeric_limits<mc::t_float>::max() << '\n';
			std::cout << "   ,lowest: " << std::numeric_limits<mc::t_float>::lowest() << '\n';
			std::cout << "   ,epsilon: " << std::numeric_limits<mc::t_float>::epsilon() << '\n';
			std::cout << "   ,round_error: " << std::numeric_limits<mc::t_float>::round_error() << '\n';


			std::cout << "***\n";

			_number_of_flights_without_boundary_collisions = 0;

			std::cin.ignore();

		}
	}
};

} // mc namespace
