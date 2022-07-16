#pragma once

#include <vector>
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\IGA.h"

class IGA_1D : public IGA
{
public:
	// Constructors
	IGA_1D();
	IGA_1D(double src, double bc, Bspline&);
	IGA_1D(IGA_1D&);

	// Destructor
	~IGA_1D();

	// Member functions
	void calcStiff();
	void calcRhs();
	void calcInteriorBasis();
	void expandSol();

	// Member variables
	Bspline bspline_x;
};