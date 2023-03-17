#ifndef H_BSPLINE
#define H_BSPLINE

#include "..\include\Utilities.h"
#include "..\include\BasisFunctions.h"
#include <iostream>
#include <vector>

class Bspline
{
public:
	Bspline() {}
	Bspline(KnotVector<double>& knewKnotVector);
	Bspline(const Bspline&);

	~Bspline() {}

	Bspline& operator=(const Bspline&);

	int findSpanOfValue(const double point);
	std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(const double point);

	KnotVector<double>& getKnotvector();
	BasisFunctions &getBasisFunctions();
	int getDegree();
	int getNumberOfBasisFunctions();

	void setKnotvector(KnotVector<double> &newKnotVector);

	void plot2D(const int resolution, std::string filename);

private:
	KnotVector<double> knotVector;
	BasisFunctions basisFunctions;
};

#include "..\src\Bspline.cpp"

#endif