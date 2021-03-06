#pragma once

#include <vector>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA.h"

class IGA_2D : public IGA
{
public:
	// Constructors
	IGA_2D();
	IGA_2D(double src, double bc, Bspline&, Bspline&);
	IGA_2D(IGA_2D&);

	// Destructor
	~IGA_2D();

	// Member functions
	void calcStiff();
	void calcRhs();
	void calcInteriorBasis();
	void LUsolver(Matrix, std::vector<double>);
	void JacobiSolver(Matrix, std::vector<double>, int iters = 10);

	// Member variables
	Bspline bspline_x;
	Bspline bspline_y;
};