
// ==================================================
// ================ 1D Poisson Problem ==============
// ==================================================

#include <iostream>
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Bspline.h"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\src\Bspline.cpp"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Matrix.h"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\src\Matrix.cpp"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\IGA.h"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\src\IGA.cpp"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\IGA_1D.h"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\src\IGA_1D.cpp"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Poisson_1D.h"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\src\Poisson_1D.cpp"

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

	// - - - - - Poisson info - - - - -
	double src = 3.0;
	double bc = 0.0;
	Poisson_1D my_poisson(src, bc, bspline_x);

	// - - - - - Compute stiffness matrix and rhs - - - - - 
	my_poisson.discretize();

	// - - - - - Solve - - - - - 
	my_poisson.solve("LU"); // "LU" or "Jacobi"
	my_poisson.expandSol();

	// - - - - - Write solution data - - - - - 
	my_poisson.plotSol("Poisson1D_sol.dat");

	return 0;
}