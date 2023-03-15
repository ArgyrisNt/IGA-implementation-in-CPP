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

	~Bspline();

	Bspline& operator=(const Bspline&);
	
	void plot(int resolution);

	int findSpanOfValue(double point);
	std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(double point);

	KnotVector<double>& getKnotvector();
	int getDegree();
	BasisFunctions getBasisFunctions();
	int getNumberOfBasisFunctions();
	Bspline &getBspline_x() { return *this; }

	void setKnotvector(KnotVector<double> &newKnotVector);

private:
	KnotVector<double> knotVector;
	BasisFunctions basisFunctions;
};

#include "..\src\Bspline.cpp"

#endif