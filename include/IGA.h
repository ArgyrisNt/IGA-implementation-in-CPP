#pragma once

#include <vector>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Matrix.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Bspline.h"

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
	void JacobiSolver(Matrix, std::vector<double>, int iters = 10);
	void LUsolver(Matrix, std::vector<double>);
	void expandSol(int);

	// Member variables
	double f;
	double bc_cond;
	Matrix stiff;
	std::vector<double> rhs;
	std::vector<double> solution;
	std::vector<double> interior_basis;	
};