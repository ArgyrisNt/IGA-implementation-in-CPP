// Solve Poisson 1D problem on a Nurbs circular line

#include <iostream>
#include "..\include\Assembler_1D.h"
#include "..\include\Poisson.h"
#include "..\include\BsplineCurve.h"

int main()
{   
	// - - - - - B-spline basis - - - - - 
	int degree = 2;
	std::vector<double> values{ 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0 };
	std::vector<double> weights{ 1.0, 0.8, 0.8, 1.0 };
	KnotVector<double> knotVector(degree, values, weights);
	Bspline bspline_x(knotVector);

	// - - - - - B-spline curve - - - - -
	std::vector<Vertex<double>> controlPoints{{0.0, 0.0}, {1.0, 1.0}, {2.0, 1.0}, {3.0, 0.0}};
	BsplineCurve curve(std::vector<Bspline>{bspline_x}, controlPoints);
	for (int i = 0; i < 1; ++i)
	{
		curve.uniformRefine_x();
	}

	// - - - - - Assempler info - - - - -
	double sourceFunction = 3.0;
	std::pair<std::string, double> west = std::make_pair("Dirichlet", 0.0);
	std::pair<std::string, double> east = std::make_pair("Dirichlet", 0.0);
	BoundCond boundaryConditions(west, east);
	Assembler_1D assembler(sourceFunction, boundaryConditions, curve);
	assembler.assemble();
	assembler.enforceBoundaryConditions("Ellimination");

	// - - - - - Poisson info - - - - -
	Poisson<Assembler_1D> poisson(assembler, Solver::Jacobi);
	poisson.solve();
	std::cout << "Solution is: " << std::endl;
	std::cout << poisson.getSolution();

	// - - - - - Plot solution - - - - -
	int resolution = 100;
	curve.plotVectorOnEntity(resolution, poisson.getSolution(), "solution.csv");

	return 0;
}