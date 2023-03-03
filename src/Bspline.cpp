#include "..\include\Bspline.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

Bspline::Bspline()
{
	knotvector = {};
	degree = 1;
	nOF = 0;
}

Bspline::Bspline(int deg, std::vector<double> &U, std::vector<double> &W)
{
	weights = W;
	knotvector = U;
	degree = deg;
	nOF = knotvector.size() - degree - 1;
	getUniqueKnots();
}

Bspline::Bspline(const Bspline &bas)
{
	knotvector = bas.knotvector;
	weights = bas.weights;
	degree = bas.degree;
	nOF = bas.nOF;
	getUniqueKnots();
}

Bspline::Bspline(const double start, const double end, int deg, int elems, std::vector<double> &W)
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

	weights = W;
	knotvector = U;
	degree = deg;
	nOF = knotvector.size() - degree - 1;
	getUniqueKnots();
}

Bspline::~Bspline() {}

Bspline &Bspline::operator=(const Bspline &old)
{
	knotvector = old.knotvector;
	weights = old.weights;
	degree = old.degree;
	nOF = old.nOF;
	getUniqueKnots();

	return *this;
}

std::pair<std::vector<double>, std::vector<double>> Bspline::calcGaussPts(int N, const double a, const double b)
{
	assert(N > 0);
	if (N > 5)
		N = 5;
	std::vector<double> GS_pts_temp;
	std::vector<double> GS_wgts;

	// Compute Gauss points in interval [0,1]
	switch (N)
	{
	case 1:
		GS_pts_temp = {0.0};
		GS_wgts = {2.0};
		break;
	case 2:
		GS_pts_temp = {-0.57735, 0.57735};
		GS_wgts = {1.0, 1.0};
		break;
	case 3:
		GS_pts_temp = {0.0, -0.774597, 0.774597};
		GS_wgts = {0.888889, 0.555556, 0.555556};
		break;
	case 4:
		GS_pts_temp = {-0.861136, -0.339981, 0.339981, 0.861136};
		GS_wgts = {0.347855, 0.652145, 0.652145, 0.347855};
		break;
	case 5:
		GS_pts_temp = {-0.90618, -0.538469, 0.0, 0.538469, 0.90618};
		GS_wgts = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927};
		break;
	default:
		std::cout << "Invalid dimension. Valid dimensions are 1,2,3,4,5." << std::endl;
		throw std::invalid_argument("Invalid dimension");
		break;
	}

	// Convert to interval [a,b]
	std::vector<double> GS_pts;
	for (int i = 0; i < GS_pts_temp.size(); i++)
	{
		GS_pts.push_back((a * (1 - GS_pts_temp[i]) + b * (1 + GS_pts_temp[i])) / 2);
	}

	return std::make_pair(GS_pts, GS_wgts);
}

int Bspline::findSpan(const double x)
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

std::pair<std::vector<double>, std::vector<double>> Bspline::eval(const double x, bool all)
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

	if (!all)
	{
		// if x in [ui,ui+1), then active functions are Ni-p,...,Ni
		std::vector<double> active_basis_func;
		std::vector<double> active_ders;
		int span = findSpan(x);
		double sum = 0.0, sum_der = 0.0;
		for (int i = span - degree; i <= span; i++)
		{
			active_basis_func.push_back(basis_func[i]);
			active_ders.push_back(ders[i]);
			sum += basis_func[i] * weights[i];
			sum_der += ders[i] * weights[i];
		}

		std::vector<double> nurbs, nurbs_ders;
		for (int i = 0; i < active_basis_func.size(); i++)
		{
			nurbs.push_back((active_basis_func[i] * weights[i]) / sum);
			nurbs_ders.push_back((active_ders[i] * weights[i] * sum - active_basis_func[i] * weights[i] * sum_der) / (sum * sum));
		}

		return std::make_pair(nurbs, nurbs_ders);
	}
	else
	{
		return std::make_pair(basis_func, ders);
	}
}

void Bspline::getUniqueKnots()
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

void Bspline::plot_basis()
{
	std::vector<std::vector<double>> all_res;
	std::vector<double> all_x;
	for (int i = 0; i <= 100; i++)
	{
		double x = i * 0.01 * knotvector[knotvector.size() - 1];
		all_x.push_back(x);
		std::vector<double> res = eval(x, true).first;
		all_res.push_back(res);
	}

	std::string filename("basis.dat");
	std::ofstream my_file(filename);
	my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
	my_file << "zone t= " << "\"1\"" << ",i=" << all_res.size() * all_res[0].size() << ",j=" << all_res.size() << "\n";

	for (int i = 0; i < all_res[0].size(); i++)
	{
		for (int j = 0; j < all_res.size(); j++)
		{
			my_file << all_x[j] << " " << all_res[j][i] << "\n";
		}
	}
	my_file.close();
}

std::vector<std::vector<double>> BsplineCurve::evaluate(bool plot)
{
	// Create B-spline curve
	std::vector<std::vector<double>> C;
	double left_limit_x = bspline_x.getKnotvector()[0];
	double right_limit_x = bspline_x.getKnotvector()[bspline_x.getKnotvector().size() - 1];

	std::string filename("curve.dat");
	std::ofstream my_file(filename);
	if (plot)
	{
		my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
   		my_file << "zone t= " << "\"1\"" << ",i=" << 101 << ",j=" << 101 << "\n";
	}

	for (int i = (int) (left_limit_x); i <= 100; i++)
	{
		double i_step = left_limit_x + (double)(i) * ((right_limit_x - left_limit_x) / 100.0);
		int span = bspline_x.findSpan(i_step);
		std::vector<double> bVal = bspline_x.eval(i_step).first;

		double coord_x = 0.0, coord_y = 0.0;
		for (int kk = 0; kk < bVal.size(); kk++)
		{
			coord_x += bVal[kk] * ctrlPts[span - bspline_x.getDegree() + kk][0];
			coord_y += bVal[kk] * ctrlPts[span - bspline_x.getDegree() + kk][1];
		}
		C.push_back({coord_x, coord_y});
		if (plot) my_file << coord_x << " " << coord_y << "\n";
	}
	my_file.close();

	return C;
}

BsplineSurface &BsplineSurface::operator=(const BsplineSurface &old)
{
	bspline_x = old.bspline_x;
	bspline_y = old.bspline_y;
	ctrlPts = old.ctrlPts;

	return *this;
}

std::vector<std::vector<double>> BsplineSurface::evaluate(bool plot)
{
	std::vector<std::vector<double>> C;
	double left_limit_x = getBspline_x().getKnotvector()[0];
	double right_limit_x = getBspline_x().getKnotvector()[getBspline_x().getKnotvector().size() - 1];
	double left_limit_y = getBspline_y().getKnotvector()[0];
	double right_limit_y = getBspline_y().getKnotvector()[getBspline_y().getKnotvector().size() - 1];
	std::string filename("surface.dat");
	std::ofstream my_file(filename);
	if (plot)
	{
		my_file << "variables= "
				<< "\"x\""
				<< ","
				<< "\"y\""
				<< "\n";
		my_file << "zone t= "
				<< "\"1\""
				<< ",i=" << 101 << ",j=" << 101 << "\n";
	}
	for (int j = (int)(left_limit_y); j <= 100; j++)
	{
		for (int i = (int)(left_limit_x); i <= 100; i++)
		{
			double i_step = left_limit_x + (double)(i) * ((right_limit_x - left_limit_x) / 100.0);
			double j_step = left_limit_y + (double)(j) * ((right_limit_y - left_limit_y) / 100.0);
			int span_i = getBspline_x().findSpan(i_step);
			std::vector<double> bVal_i = getBspline_x().eval(i_step).first;
			int span_j = getBspline_y().findSpan(j_step);
			std::vector<double> bVal_j = getBspline_y().eval(j_step).first;

			double coord_x = 0.0, coord_y = 0.0;
			for (int kkx = 0; kkx < bVal_i.size(); kkx++)
			{
				for (int kky = 0; kky < bVal_j.size(); kky++)
				{
					int i1 = span_i - getBspline_x().getDegree() + kkx;
					int i2 = span_j - getBspline_y().getDegree() + kky;
					int my = i1 * getBspline_y().getNOF() + i2;
					coord_x += bVal_i[kkx] * bVal_j[kky] * getCtrlPts()[my][0];
					coord_y += bVal_i[kkx] * bVal_j[kky] * getCtrlPts()[my][1];
				}
			}
			C.push_back({coord_x, coord_y});
			if (plot)
				my_file << coord_x << " " << coord_y << "\n";
		}
	}
	return C;
}

int Bspline::findSpanInKnotvector(std::vector<double> _knotvector, const double x)
{
	int cnt = 0;
	int index = 0;
	if (x == _knotvector[_knotvector.size() - 1])
	{
		return (_knotvector.size() - degree - 2);
	}
	for (int i = 0; i < _knotvector.size() - 1; i++)
	{
		double start = _knotvector[i];
		double end = _knotvector[i + 1];
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

void BsplineSurface::uniformRefine_x()
{
	std::vector<std::vector<std::vector<double>>> all_new_ctrlPts;
	std::vector<double> old_knotvector;
	for (int j = 0; j < bspline_y.getWeights().size(); j++)
	{
		std::vector<std::vector<double>> ctrlPts_x;
		for (int i = 0; i < ctrlPts.size(); i++)
		{
			if (i == 0 || (i % bspline_y.getWeights().size()) == 0) ctrlPts_x.push_back(ctrlPts[i + j]);
		}

		std::vector<std::vector<double>> new_ctrlPts_x;
		old_knotvector = bspline_x.getKnotvector();
		std::vector<double> new_knotvector;
		for (int i = 0; i < bspline_x.getKnotvector().size() - 1; i++)
		{
			double a = bspline_x.getKnotvector()[i];
			double b = bspline_x.getKnotvector()[i + 1];
			if (a == b) continue;
			double new_knot = (a + b) / 2.0;
			for (auto el : old_knotvector)
			{
				if (el <= new_knot) new_knotvector.push_back(el);
			}
			new_knotvector.push_back(new_knot);
			for (auto el : old_knotvector)
			{
				if (el >= new_knot) new_knotvector.push_back(el);
			}
			int span = bspline_x.findSpanInKnotvector(old_knotvector, new_knot);

			std::vector<std::vector<double>> new_ctrlPts;
			for (int ii = span - bspline_x.getDegree() + 1; ii <= span; ii++)
			{
				double a_i = (new_knot - old_knotvector[ii]) / (old_knotvector[ii + bspline_x.getDegree()] - old_knotvector[ii]);
				double coord_x = (1.0 - a_i) * ctrlPts_x[ii - 1][0] + a_i * ctrlPts_x[ii][0];
				double coord_y = (1.0 - a_i) * ctrlPts_x[ii - 1][1] + a_i * ctrlPts_x[ii][1];
				new_ctrlPts.push_back({coord_x, coord_y});
			}


			for (int ii = 0; ii <= span - bspline_x.getDegree(); ii++)
			{
				new_ctrlPts_x.push_back(ctrlPts_x[ii]);
			}
			for (int ii = 0; ii < new_ctrlPts.size(); ii++) 
			{
				new_ctrlPts_x.push_back(new_ctrlPts[ii]);
			}
			for (int ii = span; ii < ctrlPts_x.size(); ii++)
			{
				new_ctrlPts_x.push_back(ctrlPts_x[ii]);
			}
			ctrlPts_x = new_ctrlPts_x;
			new_ctrlPts_x = {};
			old_knotvector = new_knotvector;
			new_knotvector = {};
		}
		all_new_ctrlPts.push_back(ctrlPts_x);
	}

	ctrlPts = {};
	for (int i = 0; i < all_new_ctrlPts[0].size(); i++)
	{
		for (auto group : all_new_ctrlPts)
		{
			ctrlPts.push_back(group[i]);
		}
	}

	std::vector<double> new_weights(old_knotvector.size() - bspline_x.getDegree() - 1, 1.0);
	Bspline new_bspline(bspline_x.getDegree(), old_knotvector, new_weights);
	(*this).setBspline_x(new_bspline);
	(*this).setCtrlPts(ctrlPts);
}

void BsplineSurface::uniformRefine_y()
{
	std::vector<std::vector<std::vector<double>>> all_new_ctrlPts;
	std::vector<double> old_knotvector;
	for (int j = 0; j < bspline_x.getWeights().size(); j++)
	{
		std::vector<std::vector<double>> ctrlPts_x;
		for (int i = 0; i < ctrlPts.size(); i++)
		{
			if (i >= 0 && (i < bspline_y.getWeights().size()))
				ctrlPts_x.push_back(ctrlPts[i + j * bspline_y.getWeights().size()]);
		}

		std::vector<std::vector<double>> new_ctrlPts_x;
		old_knotvector = bspline_y.getKnotvector();
		std::vector<double> new_knotvector;
		for (int i = 0; i < bspline_y.getKnotvector().size() - 1; i++)
		{
			double a = bspline_y.getKnotvector()[i];
			double b = bspline_y.getKnotvector()[i + 1];
			if (a == b) continue;
			double new_knot = (a + b) / 2.0;
			for (auto el : old_knotvector)
			{
				if (el <= new_knot) new_knotvector.push_back(el);
			}
			new_knotvector.push_back(new_knot);
			for (auto el : old_knotvector)
			{
				if (el >= new_knot) new_knotvector.push_back(el);
			}
			int span = bspline_y.findSpanInKnotvector(old_knotvector, new_knot);

			std::vector<std::vector<double>> new_ctrlPts;
			for (int ii = span - bspline_y.getDegree() + 1; ii <= span; ii++)
			{
				double a_i = (new_knot - old_knotvector[ii]) / (old_knotvector[ii + bspline_y.getDegree()] - old_knotvector[ii]);
				double coord_x = (1.0 - a_i) * ctrlPts_x[ii - 1][0] + a_i * ctrlPts_x[ii][0];
				double coord_y = (1.0 - a_i) * ctrlPts_x[ii - 1][1] + a_i * ctrlPts_x[ii][1];
				new_ctrlPts.push_back({coord_x, coord_y});
			}

			for (int ii = 0; ii <= span - bspline_y.getDegree(); ii++)
			{
				new_ctrlPts_x.push_back(ctrlPts_x[ii]);
			}
			for (int ii = 0; ii < new_ctrlPts.size(); ii++)
			{
				new_ctrlPts_x.push_back(new_ctrlPts[ii]);
			}
			for (int ii = span; ii < ctrlPts_x.size(); ii++)
			{
				new_ctrlPts_x.push_back(ctrlPts_x[ii]);
			}
			ctrlPts_x = new_ctrlPts_x;
			new_ctrlPts_x = {};
			old_knotvector = new_knotvector;
			new_knotvector = {};
		}
		all_new_ctrlPts.push_back(ctrlPts_x);
	}

	ctrlPts = {};
	for (auto group : all_new_ctrlPts)
	{
		for (auto el : group) ctrlPts.push_back(el);
	}

	std::vector<double> new_weights_y(old_knotvector.size() - bspline_y.getDegree() - 1, 1.0);
	Bspline new_bspline_y(bspline_y.getDegree(), old_knotvector, new_weights_y);
	(*this).setBspline_y(new_bspline_y);
	(*this).setCtrlPts(ctrlPts);
}
