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
	std::vector<std::vector<double>> ValuesOfBasisFunctions;
	std::vector<double> evaluationPoints;
	for (int i = 0; i < resolution; ++i)
	{
		double currentStep = firstKnot + (double)(i) * ((lastKnot - firstKnot) / ((double)(resolution - 1)));
		evaluationPoints.push_back(currentStep);
		std::vector<double> values = basisFunctions->evaluateAtPoint(currentStep, true).first;
		ValuesOfBasisFunctions.push_back(values);
	}

	std::ofstream plotBspline(filename);
	plotBspline << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
	plotBspline << "zone t= " << "\"1\"" << ",i=" << resolution * ValuesOfBasisFunctions[0].size() << ",j=" << resolution << "\n";

	for (int i = 0; i < ValuesOfBasisFunctions[0].size(); ++i)
	{
		for (int j = 0; j < resolution; ++j)
		{
			plotBspline << evaluationPoints[j] << " " << ValuesOfBasisFunctions[j][i] << "\n";
		}
	}
	plotBspline.close();
}
