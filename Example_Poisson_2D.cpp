
// ==================================================
// ================ 2D Poisson Problem ==============
// ==================================================

#include <iostream>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Bspline.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\Bspline.cpp"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Matrix.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\Matrix.cpp"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\IGA.cpp"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA_2D.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\IGA_2D.cpp"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Poisson_2D.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\Poisson_2D.cpp"

// - - - - - - - - - - - - - TO DO - - - - - - - - - - - - -
//  
// • bspline_x != bspline_y has problem in plot_sol function
// • Neumann boundary conditions
// • Plot B-spline basis
// 

int main()
{
	// - - - - - B-spline basis on x-direction - - - - - 
	int p = 2;
	std::vector<double> U{ 0.0,0.0,0.0,1.0,2.0,3.0,3.0,3.0 };
	Bspline bspline_x(p, U);

	// - - - - - B-spline basis on y-direction - - - - - 
	Bspline bspline_y(bspline_x);

	// - - - - - Poisson info - - - - -
	double src = 3.0;
	double bc = 0.0;
	Poisson_2D my_poisson(src, bc, bspline_x, bspline_y);

	// - - - - - Compute stiffness matrix and rhs - - - - - 
	my_poisson.discretize();

	// - - - - - Solve - - - - - 
	my_poisson.LUsolver(my_poisson.stiff, my_poisson.rhs);

	// - - - - - Write solution data - - - - - 
	my_poisson.plotSol("Poisson2D_sol.dat", "Poisson2D_mesh.dat");

	return 0;
}