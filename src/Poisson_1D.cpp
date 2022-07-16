#include <iostream>
#include <fstream>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Poisson_1D.h"

Poisson_1D::Poisson_1D() {}

Poisson_1D::Poisson_1D(double src, double bc, Bspline& bas)
{
	f = src;
	bc_cond = bc;
	bspline_x = bas;
}

Poisson_1D::Poisson_1D(Poisson_1D& bas)
{
	f = bas.f;
	bc_cond = bas.bc_cond;
	stiff = bas.stiff;
	rhs = bas.rhs;
	solution = bas.solution;
	interior_basis = bas.interior_basis;
	bspline_x = bas.bspline_x;
}

Poisson_1D::~Poisson_1D() {}

void Poisson_1D::discretize()
{
	calcStiff();
	calcRhs();
	calcInteriorBasis();
	applyBound();
}

void Poisson_1D::applyBound()
{
	Matrix new_stiff(stiff.rows - 2, stiff.cols - 2);
	std::vector<double> new_rhs;
	int i = 0;
	for (auto it1 = interior_basis.begin(); it1 != interior_basis.end(); it1++)
	{
		int j = 0;
		for (auto it2 = interior_basis.begin(); it2 != interior_basis.end(); it2++)
		{
			new_stiff.setValue(i, j, stiff(*it1, *it2));
			j++;
		}
		new_rhs.push_back(rhs[*it1]);
		i++;
	}
	stiff = new_stiff;
	rhs = new_rhs;
}

void Poisson_1D::plotSol(std::string filename)
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
		my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "," << "\"sol\"" << "\n";
		my_file << "zone t= " << "\"1\"" << ",i=" << right_limit_x + 1 - left_limit_x << ",j=" << right_limit_x + 1 - left_limit_x << "\n";
		for (auto pt : C)
		{
			my_file << pt[0] << " " << 5.0 << " " << pt[1] << "\n";
		}
		my_file.close();
	}
}