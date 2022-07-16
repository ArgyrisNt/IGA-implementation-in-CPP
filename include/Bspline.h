#pragma once
#include <vector>

class Bspline
{
public:
	// Constructors
	Bspline(); // default constructor
	Bspline(int, std::vector<double>&);
	Bspline(Bspline&); // copy constructor
	Bspline(double start, double end, int degree, int elems);

	// Destructor
	~Bspline();

	// Member functions
	void calcGaussPts(int, double, double); // compute Gauss points and their weigths
	int findSpan(double); // compute the span that a value belongs
	std::vector<double> eval(double x); // eval all basis functions at x
	std::vector<double> ders_eval(double x); // eval all basis functions' derivatives at x

	// Member variables
	std::vector<double> knotvector;
	int degree;
	std::vector<double> GS_pts;
	std::vector<double> GS_wgts;
	std::vector<double> knots;
	int nOF;

private:
	// Member functions
	void getKnots1D(); // compute discrete knots
};