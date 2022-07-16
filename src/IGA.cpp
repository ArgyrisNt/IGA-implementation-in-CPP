#include <iostream>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Bspline.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Matrix.h"
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA.h"

IGA::IGA() {}

IGA::IGA(double src, double bc)
{
	f = src;
	bc_cond = bc;
}

IGA::IGA(IGA& bas)
{
	f = bas.f;
	bc_cond = bas.bc_cond;
	stiff = bas.stiff;
	rhs = bas.rhs;
	solution = bas.solution;
	interior_basis = bas.interior_basis;
}

IGA::~IGA() {}

void IGA::solve(std::string solver, int iters = 10)
{
	if (solver == "Jacobi")
	{
		solution = stiff.Jacobi_iterator(rhs, iters);
	}
	else if (solver == "LU")
	{
		std::vector<Matrix> LU = stiff.LU_factor();
		Matrix L = LU[0];
		Matrix U = LU[1];
		std::vector<double> y;
		y = L.forward_Euler(rhs);
		solution = U.backward_Euler(y);
	}
}
