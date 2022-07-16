#include <iostream>
#include <fstream>
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Diffusion_1D.h"

Diffusion_1D::Diffusion_1D() {}

Diffusion_1D::Diffusion_1D(double src, double bc, Bspline& bas, double k, double delta)
{
	f = src;
	bc_cond = bc;
	bspline_x = bas;
	coef = k;
	Dt = delta;
}

Diffusion_1D::Diffusion_1D(Diffusion_1D& bas)
{
	f = bas.f;
	bc_cond = bas.bc_cond;
	stiff = bas.stiff;
	rhs = bas.rhs;
	solution = bas.solution;
	interior_basis = bas.interior_basis;
	bspline_x = bas.bspline_x;
	coef = bas.coef;
	Dt = bas.Dt;
}

Diffusion_1D::~Diffusion_1D() {}

void Diffusion_1D::discretize()
{
	calcStiff();
	calcMass();
	calcRhs();
	sysMat = mass + stiff * (coef * Dt);
	for (int i = 0; i < rhs.size(); i++)
	{
		rhs[i] = rhs[i] * Dt;
	}
	calcInteriorBasis();
	applyBound();
}

void Diffusion_1D::applyBound()
{
	Matrix new_sysMatrix(sysMat.rows - 2, sysMat.cols - 2);
	std::vector<double> new_rhs;
	int i = 0;
	for (auto it1 = interior_basis.begin(); it1 != interior_basis.end(); it1++)
	{
		int j = 0;
		for (auto it2 = interior_basis.begin(); it2 != interior_basis.end(); it2++)
		{
			new_sysMatrix.setValue(i, j, sysMat(*it1, *it2));
			j++;
		}
		new_rhs.push_back(rhs[*it1]);
		i++;
	}
	sysMat = new_sysMatrix;
	rhs = new_rhs;
}

void Diffusion_1D::plotSol(std::string filename)
{
	// Create B-spline curve
	std::vector<std::vector<double>> C;
	std::vector<double> linspaced;
	double start = bspline_x.knotvector[0];
	double end = bspline_x.knotvector[bspline_x.knotvector.size() - 1];
	double delta = (end - start) / (bspline_x.nOF - 1);
	for (int i = 0; i < bspline_x.nOF - 1; i++)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end);

	std::vector<std::vector<double>> Points;
	for (int i = 0; i < linspaced.size(); i++)
	{
		Points.push_back({ linspaced[i],solution[i] });
	}

	int left_limit_x = bspline_x.knotvector[0] * 10;
	int right_limit_x = bspline_x.knotvector[bspline_x.knotvector.size() - 1] * 10;
	for (int i = left_limit_x; i <= right_limit_x; i++)
	{
		int span = bspline_x.findSpan(i * 0.1);
		std::vector<double> bVal = bspline_x.eval(i * 0.1);
		if (span + 1 > Points.size())
		{
			C.push_back(Points[Points.size() - 1]);
		}
		else
		{
			double val_x = 0.0, val_y = 0.0;
			for (int j = 0; j < bVal.size(); j++)
			{
				val_x += bVal[j] * Points[span - bspline_x.degree + j][0];
				val_y += bVal[j] * Points[span - bspline_x.degree + j][1];
			}
			C.push_back({ val_x,val_y });
		}
	}

	std::ofstream my_file(filename);
	if (my_file.is_open())
	{
		my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
		my_file << "zone t= " << "\"1\"" << ",i=" << right_limit_x + 1 - left_limit_x << ",j=" << right_limit_x + 1 - left_limit_x << "\n";
		for (auto pt : C)
		{
			my_file << pt[0] << " " << pt[1] << "\n";
		}
		my_file.close();
	}
}

void Diffusion_1D::calcMass()
{
	// Assemble mass matrix
	Matrix M(bspline_x.nOF, bspline_x.nOF, 0.0);
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

					double wvol = bspline_x.GS_wgts[g1];

					v += (bi_0 * bj_0) * wvol;
				}

				double temp = M(i1, j1) + v;
				M.setValue(i1, j1, temp);
			}
		}
	}
	mass = M;
}

std::vector<double> Diffusion_1D::nextStep()
{
	std::vector<double> b(rhs.size(), 0.0);
	std::vector<double> temp = mass * solution;
	for (int i = 0; i < rhs.size(); i++)
	{
		b[i] = rhs[i] + temp[i];
	}

	return b;
}

void Diffusion_1D::applyInitCond()
{
	std::vector<double> linspaced;
	double start = bspline_x.knotvector[0];
	double end = bspline_x.knotvector[bspline_x.knotvector.size() - 1];
	double delta = (end - start) / (bspline_x.nOF - 1);
	for (int i = 0; i < bspline_x.nOF - 1; i++)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end);
	for (int i = 0; i < linspaced.size(); i++)
	{
		if (linspaced[i] >= start && linspaced[i] <= end / 2)
		{
			solution.push_back(2.0 * linspaced[i]);
		}
		else if (linspaced[i] > end / 2 && linspaced[i] <= end)
		{
			solution.push_back(2.0 - 2.0 * linspaced[i]);
		}
	}
}
