// Solve Poisson 1D problem on a Nurbs circular line

#include <iostream>
#include "..\IGA.h"

int main()
{   
	// - - - - - B-spline basis on x-direction - - - - - 
	int p = 2;
	std::vector<double> U{ 0.0,0.0,0.0,0.5,0.5,1.0,1.0,1.0 };
	std::vector<double> W{ 1.0,sqrt(2.0)/2.0,1.0,sqrt(2.0)/2.0,1.0 };
	Bspline bspline_x(p, U, W);
	//bspline_x.plot_basis();

	// - - - - - B-spline curve - - - - -
	std::vector<std::vector<double>> ctrlPts{{2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}};
	BsplineCurve curve(bspline_x, ctrlPts);

	// - - - - - Assempler info - - - - -
	double src = 3.0;
	BoundCond _bc("Dirichlet", "Dirichlet", 0.0, 0.0);
	Assembler_1D ass1(src, _bc, curve);
	ass1.assemble();

	// - - - - - Enforce boundary conditions - - - - -
	std::string mode("Ellimination");
	ass1.enforceBoundary(mode);

	// - - - - - Poisson info - - - - -
	Poisson<Assembler_1D> poisson(ass1, Solver::ConjugateGradient);
	poisson.setSolution(poisson.getSolver()->solve(ass1.getStiff(), ass1.getRhs()));
	poisson.constructSol();
	std::cout << poisson.getSolution();

	// - - - - - Write solution data - - - - -
	poisson.plotSol("curve.dat", "sol.dat");

	return 0;
}