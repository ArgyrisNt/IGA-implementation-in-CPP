#include <iostream>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA_1D.h"

IGA_1D::IGA_1D() {}

IGA_1D::IGA_1D(double src, double bc, Bspline& bas)
{
	f = src;
	bc_cond = bc;
	bspline_x = bas;
}

IGA_1D::IGA_1D(IGA_1D& bas)
{
	f = bas.f;
	bc_cond = bas.bc_cond;
	stiff = bas.stiff;
	rhs = bas.rhs;
	solution = bas.solution;
	interior_basis = bas.interior_basis;
	bspline_x = bas.bspline_x;
}

IGA_1D::~IGA_1D() {}

void IGA_1D::calcStiff()
{
	// Assemble stiffness matrix
	Matrix A(bspline_x.nOF, bspline_x.nOF, 0.0);
	int N = bspline_x.knots.size() - 1; // Number of elements
	for (int ie1 = 0; ie1 <= N; ie1++)
	{
		int i_span_1 = bspline_x.findSpan(bspline_x.knots[ie1]);
		for (int il_1 = 0; il_1 < bspline_x.degree + 1; il_1++)
		{
			for (int jl_1 = 0; jl_1 < bspline_x.degree + 1; jl_1++)
			{
				int i1 = i_span_1 - bspline_x.degree + il_1;
				int j1 = i_span_1 - bspline_x.degree + jl_1;

				double v = 0.0;
				bspline_x.calcGaussPts(bspline_x.degree + 3, bspline_x.knots[ie1], bspline_x.knots[ie1 + 1]);
				for (int g1 = 0; g1 < bspline_x.GS_pts.size(); g1++)
				{
					std::vector<double> bVal = bspline_x.eval(bspline_x.GS_pts[g1]);
					double bi_0 = bVal[il_1];
					double bj_0 = bVal[jl_1];

					std::vector<double> gradVal = bspline_x.ders_eval(bspline_x.GS_pts[g1]);
					double bi_x = gradVal[il_1];
					double bj_x = gradVal[jl_1];

					double wvol = bspline_x.GS_wgts[g1];

					v += (bi_x * bj_x) * wvol;
				}

				double temp = A(i1, j1) + v;
				A.setValue(i1, j1, temp);
			}
		}
	}
	stiff = A;
}

void IGA_1D::calcRhs()
{
	// Assemble rhs vector
	std::vector<double> b(bspline_x.nOF, 0.0);
	int N = bspline_x.knots.size() - 1;
	for (int ie1 = 0; ie1 <= N; ie1++)
	{
		int i_span_1 = bspline_x.findSpan(bspline_x.knots[ie1]);
		for (int il_1 = 0; il_1 < bspline_x.degree + 1; il_1++)
		{
			int i1 = i_span_1 - bspline_x.degree + il_1;

			double v = 0.0;
			bspline_x.calcGaussPts(bspline_x.degree + 3, bspline_x.knots[ie1], bspline_x.knots[ie1 + 1]);
			for (int g1 = 0; g1 < bspline_x.GS_pts.size(); g1++)
			{
				std::vector<double> bVal = bspline_x.eval(bspline_x.GS_pts[g1]);
				double bi_0 = bVal[il_1];

				std::vector<double> gradVal = bspline_x.ders_eval(bspline_x.GS_pts[g1]);
				double bi_x = gradVal[il_1];

				double wvol = bspline_x.GS_wgts[g1];

				v += bi_0 * f * wvol;
			}

			b[i1] += v;
		}
	}
	rhs = b;
}

void IGA_1D::calcInteriorBasis()
{
	// Apply boundary conditions
	interior_basis = {};
	for (int i = 0; i < bspline_x.nOF; i++)
	{
		if (i == 0 || i == bspline_x.nOF - 1)
		{
			continue;
		}
		else
		{
			interior_basis.push_back(i);
		}
	}
}

void IGA_1D::LUsolver(Matrix A, std::vector<double> b)
{
	IGA::LUsolver(A, b);
	expandSol(bspline_x.nOF);
}

void IGA_1D::JacobiSolver(Matrix A, std::vector<double> b, int iters)
{
	IGA::JacobiSolver(A, b, iters);
	expandSol(bspline_x.nOF);
}