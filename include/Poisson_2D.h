#pragma once

#include <vector>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA_2D.h"

class Poisson_2D : public IGA_2D
{
public:
	// Constructors
	Poisson_2D();
	Poisson_2D(double src, double bc, Bspline&, Bspline&);
	Poisson_2D(Poisson_2D&);

	// Destructor
	~Poisson_2D();

	// Member functions
	void discretize();
	void applyBound();
	void plotSol(std::string, std::string);
};