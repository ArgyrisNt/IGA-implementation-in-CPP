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

void IGA::JacobiSolver(Matrix A, std::vector<double> b, int iters)
{
	solution = A.Jacobi_iterator(b, iters);
}

void IGA::LUsolver(Matrix A, std::vector<double> b)
{
	std::vector<Matrix> LU = A.LU_factor();
	Matrix L = LU[0];
	Matrix U = LU[1];
	std::vector<double> y;
	y = L.forward_Euler(b);
	solution = U.backward_Euler(y);
}

void IGA::expandSol(int x)
{
	std::vector<double> final_sol;
	int j = 0;
	int nOF = 0;
	nOF = x;
	for (int i = 0; i < nOF; i++)
	{
		if (std::count(interior_basis.begin(), interior_basis.end(), i))
		{
			final_sol.push_back(solution[j]);
			j++;
		}
		else
		{
			final_sol.push_back(bc_cond);
		}
	}
	solution = final_sol;
}