#include "..\include\Bspline.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

Bspline::Bspline(const Bspline &bspline)
{
	knotVector = bspline.knotVector;
	basisFunctions = bspline.basisFunctions;
}

Bspline &Bspline::operator=(const Bspline &bspline)
{
	knotVector = bspline.knotVector;
	basisFunctions = bspline.basisFunctions;

	return *this;
}


void Bspline::setKnotvector(const KnotVector<double> &newKnotVector)
{
	knotVector = newKnotVector;
	basisFunctions = std::make_shared<BasisFunctions>(BasisFunctions(knotVector));
}


void Bspline::plot(const int resolution, const std::string &filename)
{
	double firstKnot = knotVector(0);
	double lastKnot = knotVector(knotVector.getSize()-1);
	std::ofstream plotBspline(filename);
	plotBspline << "X";
	for (int i = 0; i < getNumberOfBasisFunctions(); ++i)
	{
		plotBspline << "," << i;
	}
	plotBspline << "\n";
	for (int i = 0; i < resolution; ++i)
	{
		double currentStep = firstKnot + (double)(i) * ((lastKnot - firstKnot) / ((double)(resolution - 1)));
		std::vector<double> values = basisFunctions->evaluateAtPoint(currentStep, true).first;

		plotBspline << currentStep;
		for (auto it = values.begin(); it != values.end(); ++it)
		{
			plotBspline << "," << *it;
		}
		plotBspline << "\n";
	}
	plotBspline.close();
}
