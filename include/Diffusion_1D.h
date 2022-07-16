#pragma once

#include <vector>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA_1D.h"

class Diffusion_1D : public IGA_1D
{
public:
	// Constructors
	Diffusion_1D();
	Diffusion_1D(double src, double bc, Bspline&, double k, double delta);
	Diffusion_1D(Diffusion_1D&);

	// Destructor
	~Diffusion_1D();

	// Member functions
	void discretize();
	void applyBound();
	void plotSol(std::string);
	void calcMass();
	std::vector<double> nextStep();
	void applyInitCond();

	// Member variables
	Matrix mass;
	Matrix sysMat;
	double coef;
	double Dt;
};