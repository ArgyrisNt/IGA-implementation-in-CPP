#pragma once

#include <vector>
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\IGA_1D.h"

class Poisson_1D : public IGA_1D
{
public:
	// Constructors
	Poisson_1D();
	Poisson_1D(double src, double bc, Bspline&);
	Poisson_1D(Poisson_1D&);

	// Destructor
	~Poisson_1D();

	// Member functions
	void discretize();
	void applyBound();
	void plotSol(std::string);
};