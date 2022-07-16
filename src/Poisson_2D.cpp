#include <iostream>
#include <fstream>
#include "C:\Users\argir\IGA-implementation-in-CPP\include\Poisson_2D.h"

Poisson_2D::Poisson_2D() {}

Poisson_2D::Poisson_2D(double src, double bc, Bspline& bas_x, Bspline& bas_y)
{
	f = src;
	bc_cond = bc;
	bspline_x = bas_x;
	bspline_y = bas_y;
}

Poisson_2D::Poisson_2D(Poisson_2D& bas)
{
	f = bas.f;
	bc_cond = bas.bc_cond;
	bspline_x = bas.bspline_x;
	bspline_y = bas.bspline_y;
}

Poisson_2D::~Poisson_2D() {}

void Poisson_2D::discretize()
{
	calcStiff();
	calcRhs();
	calcInteriorBasis();
	applyBound();
}

void Poisson_2D::applyBound()
{
	int new_dim = interior_basis.size();
	Matrix new_stiff(new_dim, new_dim);
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

void Poisson_2D::plotSol(std::string filename1, std::string filename2)
{
	// Create B-spline surface
	std::vector<std::vector<double>> C;
	std::vector<double> linspaced_x;
	double start = bspline_x.knotvector[0];
	double end = bspline_x.knotvector[bspline_x.knotvector.size() - 1];
	double delta = (end - start) / (bspline_x.nOF - 1);
	for (int i = 0; i < bspline_x.nOF - 1; i++)
	{
		linspaced_x.push_back(start + delta * i);
	}
	linspaced_x.push_back(end);

	std::vector<double> linspaced_y;
	start = bspline_y.knotvector[0];
	end = bspline_y.knotvector[bspline_y.knotvector.size() - 1];
	delta = (end - start) / (bspline_y.nOF - 1);
	for (int i = 0; i < bspline_y.nOF - 1; i++)
	{
		linspaced_y.push_back(start + delta * i);
	}
	linspaced_y.push_back(end);

	std::vector<std::vector<double>> Points;
	int ic = 0;
	for (int i = 0; i < linspaced_x.size(); i++)
	{
		for (int j = 0; j < linspaced_y.size(); j++)
		{
			Points.push_back({ linspaced_x[i],linspaced_y[j],solution[ic] });
			ic++;
		}
	}

	int left_limit_x = bspline_x.knotvector[0] * 10;
	int right_limit_x = bspline_x.knotvector[bspline_x.knotvector.size() - 1] * 10;
	int left_limit_y = bspline_y.knotvector[0] * 10;
	int right_limit_y = bspline_y.knotvector[bspline_y.knotvector.size() - 1] * 10;
	for (int j = left_limit_y; j <= right_limit_y; j++)
	{
		for (int i = left_limit_x; i <= right_limit_x; i++)
		{
			int span_i = bspline_x.findSpan(i * 0.1);
			std::vector<double> bVal_i = bspline_x.eval(i * 0.1);
			int span_j = bspline_y.findSpan(j * 0.1);
			std::vector<double> bVal_j = bspline_y.eval(j * 0.1);
			double val_x = 0.0, val_y = 0.0, val_z = 0.0;
			for (int k = 0; k < bVal_i.size(); k++)
			{
				for (int kk = 0; kk < bVal_j.size(); kk++)
				{
					int i1 = span_i - bspline_x.degree + k;
					int i2 = span_j - bspline_y.degree + kk;
					int my = i1 * bspline_y.nOF + i2;
					val_x += bVal_i[k] * bVal_j[kk] * Points[my][0];
					val_y += bVal_i[k] * bVal_j[kk] * Points[my][1];
					val_z += bVal_i[k] * bVal_j[kk] * Points[my][2];
				}
			}
			C.push_back({ val_x,val_y,val_z });
		}
	}

	std::ofstream my_file1(filename1);
	if (my_file1.is_open())
	{
		my_file1 << "variables= " << "\"x\"" << "," << "\"y\"" << "," << "\"sol\"" << "\n";
		my_file1 << "zone t= " << "\"1\"" << ",i=" << right_limit_x + 1 - left_limit_x << ",j=" << right_limit_y + 1 - left_limit_y << "\n";
		for (auto pt : C)
		{
			my_file1 << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
		}
		my_file1.close();
	}

	std::ofstream my_file2(filename2);
	if (my_file2.is_open())
	{
		my_file2 << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
		my_file2 << "zone t= " << "\"1\"" << ",i=" << bspline_x.knots.size() << ",j=" << bspline_y.knots.size() << "\n";
		for (int j = 0; j < bspline_y.knots.size(); j++)
		{
			for (int i = 0; i < bspline_x.knots.size(); i++)
			{
				my_file2 << bspline_x.knots[i] << " " << bspline_y.knots[j] << "\n";
			}
		}
		my_file2.close();
	}
}