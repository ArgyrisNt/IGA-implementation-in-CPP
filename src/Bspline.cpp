#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Bspline.h"
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Matrix.h"
#include <algorithm>
#include <cassert>

Bspline::Bspline()
{
	knotvector = {};
	degree = 1;
	nOF = 0;
}

Bspline::Bspline(int deg, std::vector<double>& U)
{
	knotvector = U;
	degree = deg;
	nOF = knotvector.size() - degree - 1;
	getKnots1D();
}

Bspline::Bspline(Bspline& bas)
{
	knotvector = bas.knotvector;
	degree = bas.degree;
	GS_pts = bas.GS_pts;
	GS_wgts = bas.GS_wgts;
	nOF = bas.nOF;
	getKnots1D();
}

Bspline::Bspline(double start, double end, int deg, int elems)
{
	std::vector<double> U;
	for (int i = 0; i < deg; i++)
	{
		U.push_back(start);
	}
	double delta = end / elems;
	for (int i = 0; i <= elems; i++)
	{
		U.push_back(i * delta);
	}
	for (int i = 0; i < deg; i++)
	{
		U.push_back(end);
	}

	knotvector = U;
	degree = deg;
	nOF = knotvector.size() - degree - 1;
	getKnots1D();
}

Bspline::~Bspline() {}

void Bspline::calcGaussPts(int N, double a, double b)
{
	assert(N > 0);
	if (N > 5) N = 5;
	std::vector<double> GS_pts_temp;
	std::vector<double> GS_wgts_temp;

	// Compute Gauss points in interval [0,1]
	switch (N)
	{
	case 1:
		GS_pts_temp = { 0.0 };
		GS_wgts_temp = { 2.0 };
		break;
	case 2:
		GS_pts_temp = { -0.57735, 0.57735 };
		GS_wgts_temp = { 1.0, 1.0 };
		break;
	case 3:
		GS_pts_temp = { 0.0, -0.774597, 0.774597 };
		GS_wgts_temp = { 0.888889, 0.555556, 0.555556 };
		break;
	case 4:
		GS_pts_temp = { -0.861136, -0.339981, 0.339981, 0.861136 };
		GS_wgts_temp = { 0.347855, 0.652145, 0.652145, 0.347855 };
		break;
	case 5:
		GS_pts_temp = { -0.90618, -0.538469, 0.0, 0.538469, 0.90618 };
		GS_wgts_temp = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };
		break;
	default:
		std::cout << "Invalid dimension. Valid dimensions are 1,2,3,4,5." << std::endl;
		throw std::invalid_argument("Invalid dimension");
		break;
	}

	// Convert to interval [a,b]
	GS_pts = {};
	for (int i = 0; i < GS_pts_temp.size(); i++)
	{
		GS_pts.push_back((a * (1 - GS_pts_temp[i]) + b * (1 + GS_pts_temp[i])) / 2);
	}
	GS_wgts = {};
	GS_wgts = GS_wgts_temp;
}

int Bspline::findSpan(double x)
{
	int cnt = 0;
	int index = 0;
	if (x == knotvector[knotvector.size() - 1])
	{
		return (knotvector.size() - degree - 2);
	}
	for (int i = 0; i < knotvector.size() - 1; i++)
	{
		double start = knotvector[i];
		double end = knotvector[i + 1];
		if (x >= start && x < end)
		{
			index = i;
			cnt += 1;
			break;
		}
	}
	if (cnt == 0)
	{
		std::cout << "Error: Value " << x << " is not in knotvector." << std::endl;
		throw std::invalid_argument("Value does not appear in knotvector");
	}

	return index;
}

std::vector<double> Bspline::eval(double x)
{
	std::vector<double> basis_func(knotvector.size() - degree - 1);
	for (int j = 0; j < knotvector.size() - degree - 1; j++)
	{
		if (x >= knotvector[j] && x < knotvector[j + 1])
		{
			basis_func[j] = 1.0;
		}
		else
		{
			basis_func[j] = 0.0;
		}
	}

	if (x == knotvector[knotvector.size() - 1])
	{
		basis_func[basis_func.size() - 1] = 1.0;
	}

	for (int i = 1; i < degree + 1; i++)
	{
		std::vector<double> basis_func_old = basis_func;
		for (int j = 0; j < knotvector.size() - degree - 1; j++)
		{
			double nomin = x - knotvector[j];
			double denomin = knotvector[j + i] - knotvector[j];
			double coef1 = 0.0, coef2 = 0.0;
			if (denomin != 0.0)
			{
				coef1 = nomin / denomin;
			}
			nomin = knotvector[j + i + 1] - x;
			denomin = knotvector[j + i + 1] - knotvector[j + 1];
			if (denomin != 0.0)
			{
				coef2 = nomin / denomin;
			}
			if (j == knotvector.size() - degree - 2)
			{
				basis_func[j] = coef1 * basis_func_old[j];
			}
			else
			{
				basis_func[j] = coef1 * basis_func_old[j] + coef2 * basis_func_old[j + 1];
			}
		}
	}

	// if x in [ui,ui+1), then active functions are Ni-p,...,Ni
	std::vector<double> active_basis_func;
	int span = findSpan(x);
	for (int i = span - degree; i <= span; i++)
	{
		active_basis_func.push_back(basis_func[i]);
	}

	return active_basis_func;
}

std::vector<double> Bspline::ders_eval(double x)
{
	std::vector<double> basis_func(knotvector.size() - degree - 1);
	std::vector<double> ders(knotvector.size() - degree - 1);
	for (int j = 0; j < knotvector.size() - degree - 1; j++)
	{
		if (x >= knotvector[j] && x < knotvector[j + 1])
		{
			basis_func[j] = 1.0;
		}
		else
		{
			basis_func[j] = 0.0;
		}
	}

	if (x == knotvector[knotvector.size() - 1])
	{
		basis_func[basis_func.size() - 1] = 1.0;
	}

	for (int i = 1; i < degree + 1; i++)
	{
		std::vector<double> basis_func_old = basis_func;
		for (int j = 0; j < knotvector.size() - degree - 1; j++)
		{
			double nomin = x - knotvector[j];
			double denomin = knotvector[j + i] - knotvector[j];
			double coef1 = 0.0, coef2 = 0.0;
			if (denomin != 0.0)
			{
				coef1 = nomin / denomin;
			}
			nomin = knotvector[j + i + 1] - x;
			denomin = knotvector[j + i + 1] - knotvector[j + 1];
			if (denomin != 0.0)
			{
				coef2 = nomin / denomin;
			}
			if (j == knotvector.size() - degree - 2)
			{
				basis_func[j] = coef1 * basis_func_old[j];
			}
			else
			{
				basis_func[j] = coef1 * basis_func_old[j] + coef2 * basis_func_old[j + 1];
			}
			if (i == degree)
			{
				double coef3 = 0.0, coef4 = 0.0;
				if (knotvector[j + i] - knotvector[j] != 0)
				{
					coef3 = basis_func_old[j] / (knotvector[j + i] - knotvector[j]);
				}
				if (knotvector[j + i + 1] - knotvector[j + 1] != 0)
				{
					coef4 = basis_func_old[j + 1] / (knotvector[j + i + 1] - knotvector[j + 1]);
				}
				ders[j] = degree * (coef3 - coef4);
			}
		}
	}

	// if x in [ui,ui+1), then active functions are Ni-p,...,Ni
	std::vector<double> active_ders;
	int span = findSpan(x);
	for (int i = span - degree; i <= span; i++)
	{
		active_ders.push_back(ders[i]);
	}

	return active_ders;
}

void Bspline::getKnots1D()
{
	knots = {};
	for (int i = degree; i < knotvector.size() - degree; i++)
	{
		if (!std::count(knots.begin(), knots.end(), knotvector[i]))
		{
			knots.push_back(knotvector[i]);
		}
	}
}