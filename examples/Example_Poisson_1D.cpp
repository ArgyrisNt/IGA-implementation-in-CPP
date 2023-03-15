// Solve Poisson 1D problem on a Nurbs circular line

#include <iostream>
#include "..\include\Assembler_1D.h"
#include "..\include\Poisson.h"
#include "..\include\BsplineCurve.h"

int main()
{   
	// - - - - - B-spline basis on x-direction - - - - - 
	int degree = 2;
	std::vector<double> values{ 0.0,0.0,0.0,0.5,1.0,1.0,1.0 };
	std::vector<double> weights{ 1.0, 1.0, 1.0, 1.0 };
	KnotVector<double> knotVector(degree, values, weights);
	Bspline bspline_x(knotVector);
	//int resolution = 100;
	//bspline_x.plot(resolution);

	// - - - - - B-spline curve - - - - -
	std::vector<std::vector<double>> controlPoints{{0.0, 0.0}, {1.0, 1.0}, {2.0, 1.0}, {3.0, 0.0}};
	BsplineCurve curve(bspline_x, controlPoints);

	// - - - - - Assempler info - - - - -
	double sourceFunction = 3.0;
	BoundCond boundaryConditions("Dirichlet", "Dirichlet", 0.0, 0.0);
	Assembler_1D assembler(sourceFunction, boundaryConditions, curve);
	assembler.assemble();

	// - - - - - Enforce boundary conditions - - - - -
	std::string mode("Ellimination");
	assembler.enforceBoundaryConditions(mode);

	// - - - - - Poisson info - - - - -
	Poisson<Assembler_1D> poisson(assembler, Solver::Jacobi);
	poisson.setSolution(poisson.getSolver()->solve());
	std::cout << "Solution is: " << std::endl;
	std::cout << poisson.getSolution();

	// - - - - - Write solution data - - - - -
	int resolution = 100;
	curve.plot3D(resolution, poisson.getSolution(), "solution.dat");

	return 0;
}