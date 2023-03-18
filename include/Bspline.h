#ifndef H_BSPLINE
#define H_BSPLINE

#include "..\include\BasisFunctions.h"
#include <iostream>
#include <vector>

class Bspline
{
public:
	Bspline() {}
	Bspline(const KnotVector<double> &knewKnotVector);
	Bspline(const Bspline&);

	~Bspline() {}

	Bspline& operator=(const Bspline&);

	int findSpanOfValue(const double point);
	std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(const double point);

	KnotVector<double>& getKnotvector();
	BasisFunctions &getBasisFunctions();
	int getDegree();
	int getNumberOfBasisFunctions();

	void setKnotvector(const KnotVector<double> &newKnotVector);

	void plot2D(const int resolution, const std::string &filename);

private:
	KnotVector<double> knotVector;
	BasisFunctions basisFunctions;
};

#include "..\src\Bspline.cpp"

#endif