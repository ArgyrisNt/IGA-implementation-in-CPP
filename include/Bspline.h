#ifndef H_BSPLINE
#define H_BSPLINE

#include "BasisFunctions.h"
#include <iostream>
#include <vector>

class Bspline
{
public:
	Bspline() {}
	Bspline(const KnotVector<double> &newKnotVector)
		: knotVector(newKnotVector), basisFunctions(std::make_shared<BasisFunctions>(BasisFunctions(knotVector))) {}
	Bspline(const Bspline&);

	~Bspline() {}

	Bspline& operator=(const Bspline&);

	int findSpanOfValue(const double point) 
	{ return knotVector.findSpanOfValue(point); }

	std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(const double point)
	{ return basisFunctions->evaluateAtPoint(point); }

	KnotVector<double> &getKnotvector() 
	{ return knotVector; }

	const int getDegree() const 
	{ return knotVector.getDegree(); }

	const int getNumberOfBasisFunctions() const 
	{ return knotVector.getSize() - knotVector.getDegree() - 1; }

	void setKnotvector(const KnotVector<double> &newKnotVector);

	void plot(const int resolution, const std::string &filename);

private:
	KnotVector<double> knotVector;
	std::shared_ptr<BasisFunctions> basisFunctions;
};

#include "..\src\Bspline.cpp"

#endif