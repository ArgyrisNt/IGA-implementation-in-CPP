#include <iostream>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\IGA_2D.h"

IGA_2D::IGA_2D() {}

IGA_2D::IGA_2D(double src, double bc, Bspline& bas_x, Bspline& bas_y)
{
	f = src;
	bc_cond = bc;
	bspline_x = bas_x;
	bspline_y = bas_y;
}

IGA_2D::IGA_2D(IGA_2D& bas)
{
	f = bas.f;
	bc_cond = bas.bc_cond;
	stiff = bas.stiff;
	rhs = bas.rhs;
	solution = bas.solution;
	interior_basis = bas.interior_basis;
	bspline_x = bas.bspline_x;
	bspline_y = bas.bspline_y;
}

IGA_2D::~IGA_2D() {}

void IGA_2D::calcStiff()
{
	// Assemble stiffness matrix
	int Nx = bspline_x.knots.size() - 1; // Number of elements on x-direction
	int Ny = bspline_y.knots.size() - 1; // NUmber of elements on y-direction
	int nOF = bspline_x.nOF * bspline_y.nOF; // Number of total basis function
	Matrix A(nOF, nOF, 0.0);
	for (int ie1 = 0; ie1 <= Nx; ie1++)
	{
		int i_span_1 = bspline_x.findSpan(bspline_x.knots[ie1]);
		for (int ie2 = 0; ie2 <= Ny; ie2++)
		{
			int i_span_2 = bspline_y.findSpan(bspline_y.knots[ie2]);
			for (int il_1 = 0; il_1 < bspline_x.degree + 1; il_1++)
			{
				for (int il_2 = 0; il_2 < bspline_y.degree + 1; il_2++)
				{
					for (int jl_1 = 0; jl_1 < bspline_x.degree + 1; jl_1++)
					{
						for (int jl_2 = 0; jl_2 < bspline_y.degree + 1; jl_2++)
						{
							int i1 = i_span_1 - bspline_x.degree + il_1;
							int j1 = i_span_1 - bspline_x.degree + jl_1;

							int i2 = i_span_2 - bspline_y.degree + il_2;
							int j2 = i_span_2 - bspline_y.degree + jl_2;

							double v = 0.0;
							bspline_x.calcGaussPts(bspline_x.degree + 3, bspline_x.knots[ie1], bspline_x.knots[ie1 + 1]);
							std::vector<double> GS_1 = bspline_x.GS_pts;
							std::vector<double> GS_w_1 = bspline_x.GS_wgts;
							bspline_y.calcGaussPts(bspline_y.degree + 3, bspline_y.knots[ie2], bspline_y.knots[ie2 + 1]);
							std::vector<double> GS_2 = bspline_y.GS_pts;
							std::vector<double> GS_w_2 = bspline_y.GS_wgts;
							for (int g1 = 0; g1 < GS_1.size(); g1++)
							{
								for (int g2 = 0; g2 < GS_2.size(); g2++)
								{
									std::vector<double> bVal_1 = bspline_x.eval(GS_1[g1]);
									std::vector<double> bVal_2 = bspline_y.eval(GS_2[g2]);
									double bi_0 = bVal_1[il_1] * bVal_2[il_2];
									double bj_0 = bVal_1[jl_1] * bVal_2[jl_2];

									std::vector<double> gradVal_1 = bspline_x.ders_eval(GS_1[g1]);
									std::vector<double> gradVal_2 = bspline_y.ders_eval(GS_2[g2]);
									double bi_x = gradVal_1[il_1] * bVal_2[il_2];
									double bi_y = bVal_1[il_1] * gradVal_2[il_2];
									double bj_x = gradVal_1[jl_1] * bVal_2[jl_2];
									double bj_y = bVal_1[jl_1] * gradVal_2[jl_2];

									double wvol = GS_w_1[g1] * GS_w_2[g2];

									v += (bi_x * bj_x + bi_y * bj_y) * wvol;
								}
							}

							int index1 = i1 * bspline_y.nOF + i2;
							int index2 = j1 * bspline_y.nOF + j2;
							double temp = A(index1, index2) + v;
							A.setValue(index1, index2, temp);
						}
					}
				}
			}
		}
	}
	stiff = A;
}

void IGA_2D::calcRhs()
{
	// Assemble rhs vector
	int nOF = bspline_x.nOF * bspline_y.nOF;
	std::vector<double> b(nOF, 0.0);
	int Nx = bspline_x.knots.size() - 1;
	int Ny = bspline_y.knots.size() - 1;
	for (int ie1 = 0; ie1 <= Nx; ie1++)
	{
		int i_span_1 = bspline_x.findSpan(bspline_x.knots[ie1]);
		for (int ie2 = 0; ie2 <= Ny; ie2++)
		{
			int i_span_2 = bspline_y.findSpan(bspline_y.knots[ie2]);

			for (int il_1 = 0; il_1 < bspline_x.degree + 1; il_1++)
			{
				for (int il_2 = 0; il_2 < bspline_y.degree + 1; il_2++)
				{
					int i1 = i_span_1 - bspline_x.degree + il_1;
					int i2 = i_span_2 - bspline_y.degree + il_2;

					double v = 0.0;
					bspline_x.calcGaussPts(bspline_x.degree + 3, bspline_x.knots[ie1], bspline_x.knots[ie1 + 1]);
					std::vector<double> GS_1 = bspline_x.GS_pts;
					std::vector<double> GS_w_1 = bspline_x.GS_wgts;
					bspline_y.calcGaussPts(bspline_y.degree + 3, bspline_y.knots[ie2], bspline_y.knots[ie2 + 1]);
					std::vector<double> GS_2 = bspline_y.GS_pts;
					std::vector<double> GS_w_2 = bspline_y.GS_wgts;
					for (int g1 = 0; g1 < GS_1.size(); g1++)
					{
						for (int g2 = 0; g2 < GS_2.size(); g2++)
						{
							std::vector<double> bVal_1 = bspline_x.eval(GS_1[g1]);
							std::vector<double> bVal_2 = bspline_y.eval(GS_2[g2]);
							double bi_0 = bVal_1[il_1] * bVal_2[il_2];

							double wvol = GS_w_1[g1] * GS_w_2[g2];

							v += bi_0 * f * wvol;
						}
					}

					int index = i1 * bspline_y.nOF + i2;
					b[index] += v;
				}
			}
		}
	}
	rhs = b;
}

void IGA_2D::calcInteriorBasis()
{
	// Apply boundary conditions
	interior_basis = {};
	int nOF = bspline_x.nOF * bspline_y.nOF;
	for (int i = 0; i < nOF; i++)
	{
		if (i >= 0 && i <= bspline_x.nOF - 1)
		{
			continue;
		}
		else if (i % bspline_x.nOF == 0 || i % bspline_x.nOF == bspline_x.nOF - 1)
		{
			continue;
		}
		else if (i >= nOF - bspline_x.nOF)
		{
			continue;
		}
		else
		{
			interior_basis.push_back(i);
		}
	}
}

void IGA_2D::LUsolver(Matrix A, std::vector<double> b)
{
	IGA::LUsolver(A, b);
	expandSol(bspline_x.nOF * bspline_y.nOF);
}

void IGA_2D::JacobiSolver(Matrix A, std::vector<double> b, int iters)
{
	IGA::JacobiSolver(A, b, iters);
	expandSol(bspline_x.nOF * bspline_y.nOF);
}