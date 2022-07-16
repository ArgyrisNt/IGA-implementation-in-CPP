#pragma once

#include <vector>
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Matrix.h"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Bspline.h"

class IGA
{
public:
	// Constructors
	IGA(); // Default constructor
	IGA(double src, double bc);
	IGA(IGA&); // Copy constructor

	// Destructor
	~IGA();

	// Member functions
	void solve(std::string, int); // Solve the linear problem Ax=b

	// Member variables
	double f;
	double bc_cond;
	Matrix stiff;
	std::vector<double> rhs;
	std::vector<double> solution;
	std::vector<double> interior_basis;	
};