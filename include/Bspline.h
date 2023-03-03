#pragma once

#include <iostream>
#include <vector>

class Bspline
{
public:
	// Constructors
	Bspline(); // default constructor
	Bspline(int, std::vector<double>&, std::vector<double>&); // initialize a basis by defining a degree and a knotvector
	Bspline(const Bspline&); // copy constructor
	Bspline(const double start, const double end, int degree, int elems, std::vector<double>& W); // initialize a basis by defining the start and the end of the knotvector, the degree and the number of elements

	// Destructor
	~Bspline();

	// Overload operator
	Bspline& operator=(const Bspline&);

	// Member functions
	std::pair<std::vector<double>, std::vector<double>> calcGaussPts(int, const double, const double); // compute Gauss points and their weigths
	int findSpan(const double); // compute the span that a value belongs
	std::pair<std::vector<double>, std::vector<double>> eval(const double, bool all = false); // eval all active basis functions at x (all basis functions if all=true)
	void plot_basis();
	int findSpanInKnotvector(std::vector<double> _knotvector, const double x);

	// Member setter functions
	void setWeights(std::vector<double> _weights) 
	{ 
		bool invalid = false;
		for (auto el : _weights)
		{
			if (el < 0 || el > 1)
			{
				invalid = true;
				break;
			}
		}
		if (!invalid) weights = _weights; 
		else
		{
			std::cout << "Invalid weights. Valid weights are within [0,1]." << std::endl;
		    throw std::invalid_argument("Invalid weights");
		}
	}
	void setKnotvector(std::vector<double> new_knotvector) { knotvector = new_knotvector; }

	// Member getter functions
	std::vector<double>& getKnotvector() { return knotvector; };
	const int getDegree() const { return degree; }
	const int getNOF() const { return nOF; }
	std::vector<double>& getWeights() { return weights; }

	// Member variables
	std::vector<double> knots;

private:
	// Member local functions
	void getUniqueKnots(); // compute discrete knots

	// Member variables
	std::vector<double> knotvector;
	int degree;
	int nOF;
	std::vector<double> weights;
};

class BsplineCurve
{
public:
	BsplineCurve(Bspline _bspline_x, std::vector<std::vector<double>> _ctrlPts)
	{
		bspline_x = _bspline_x;
		ctrlPts = _ctrlPts;
	}

	~BsplineCurve() {}

	std::vector<std::vector<double>> evaluate(bool plot = true);

	Bspline bspline_x;
	std::vector<std::vector<double>> ctrlPts;
};

class BsplineSurface
{
public:
	BsplineSurface(const Bspline &_bspline_x, const Bspline &_bspline_y, std::vector<std::vector<double>> _ctrlPts)
	{
		bspline_x = _bspline_x;
		bspline_y = _bspline_y;
		ctrlPts = _ctrlPts;
	}

	~BsplineSurface() {}

	BsplineSurface &operator=(const BsplineSurface &);

	Bspline &getBspline_x() { return bspline_x; }
	Bspline &getBspline_y() { return bspline_y; }
	std::vector<std::vector<double>> &getCtrlPts() { return ctrlPts; };

	void setCtrlPts(std::vector<std::vector<double>> &_ctrlPts) { ctrlPts = _ctrlPts; }
	void setBspline_x(Bspline &_bspline_x) { bspline_x = _bspline_x; }
	void setBspline_y(Bspline &_bspline_y) { bspline_y = _bspline_y; }

	std::vector<std::vector<double>> evaluate(bool plot = true);
	void uniformRefine_x();
	void uniformRefine_y();

private:
	Bspline bspline_x;
	Bspline bspline_y;
	std::vector<std::vector<double>> ctrlPts;
};