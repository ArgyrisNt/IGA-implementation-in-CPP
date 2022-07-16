
// ==================================================
// ============== 1D Diffusion Problem ==============
// ==================================================

#include <iostream>
#include <string>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Bspline.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\Bspline.cpp"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Matrix.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\Matrix.cpp"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\IGA.cpp"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA_1D.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\IGA_1D.cpp"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Diffusion_1D.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\src\Diffusion_1D.cpp"

// - - - - - - - - - TO DO - - - - - - - -
// 
// • Initial condition should be an input
// • Error computation
// 

int main()
{
	// - - - - - B-spline basis on x-direction - - - - - 
	double start = 0.0;
	double end = 1.0;
	int p = 3;
	int numElems = 200;
	Bspline bspline_x(start, end, p, numElems);

	// - - - - - Diffusion info - - - - -
	double src = 0.0;
	double bc = 0.0;
	double coef = 1.0;
	double t_end = 0.05; // 1
	double numSteps = 10.0; // 200
	double Dt = t_end / numSteps;
	Diffusion_1D my_diffusion(src, bc, bspline_x, coef, Dt);

	// - - - - - Compute stiffness matrix, mass matrix and rhs - - - - - 
	my_diffusion.discretize();

	// - - - - - Apply and plot initial condition - - - - -
	my_diffusion.applyInitCond();
	my_diffusion.plotSol("0sol.dat");

	// - - - - - Solve - - - - - 
	for (int t = 0; t < numSteps; t++) // until 6
	{
		std::cout << std::endl << "---------------- " << t + 1 << " step ----------------";
		std::vector<double> b = my_diffusion.nextStep(); // build next rhs
		my_diffusion.solution = my_diffusion.sysMat.Jacobi_iterator(b, 100);
		my_diffusion.expandSol();
		my_diffusion.plotSol(std::to_string(t + 1) + "sol.dat");
	}

	return 0;
}