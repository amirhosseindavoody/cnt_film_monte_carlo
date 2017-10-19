#include <iostream>
#include <cmath>

#include "energy_dispersion.h"

#include "utility.h"

namespace mc
{

// constructor
energy_dispersion::energy_dispersion()
{
	t_float a = 5.431e-10; // silicon lattice constant as a typical value
	t_float max_k = 2.*pi/a;

	_n1 = 101;
	_n2 = 101;
	_n3 = 101; 

	t_float dk = max_k/t_float(_n1-1);

	_energy = mc::v3d(_n1,mc::v2d(_n2,mc::v1d(_n3)));
	_k1 = mc::v3d(_n1,mc::v2d(_n2,mc::v1d(_n3)));
	_k2 = mc::v3d(_n1,mc::v2d(_n2,mc::v1d(_n3)));
	_k3 = mc::v3d(_n1,mc::v2d(_n2,mc::v1d(_n3)));

	for (int i=0; i<(_n1-1)/2; ++i)
	{
		for (int j=0; j<(_n2-1)/2; ++j)	
		{
			for (int k=0; k<(_n3-1)/2; ++k)
			{
				t_float p2 = std::pow(hbar*dk,2)*t_float(i*i+j*j+k*k);
				_energy[i][j][k] = p2/(2.*mc::elec_mass);
				_k1[i][j][k] = t_float(i)*dk;
				_k2[i][j][k] = t_float(j)*dk;
				_k3[i][j][k] = t_float(k)*dk;
			}
		}
	}
};

} // mc namespace