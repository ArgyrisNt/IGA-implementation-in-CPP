#include "..\include\Bspline.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

Bspline::Bspline(KnotVector<double> &newKnotVector)
{
	knotVector = newKnotVector;
	basisFunctions.setKnotVector(knotVector);
}

Bspline::Bspline(const Bspline &bspline)
{
	knotVector = bspline.knotVector;
	basisFunctions = bspline.basisFunctions;
}

Bspline::~Bspline() {}

Bspline &Bspline::operator=(const Bspline &bspline)
{
	knotVector = bspline.knotVector;
	basisFunctions = bspline.basisFunctions;

	return *this;
}



KnotVector<double> &Bspline::getKnotvector()
{
	return knotVector;
}

int Bspline::getDegree()
{
	return knotVector.getDegree();
}

BasisFunctions Bspline::getBasisFunctions()
{
	return basisFunctions;
}

int Bspline::getNumberOfBasisFunctions()
{
	return getKnotvector().getSize() - getKnotvector().getDegree() - 1;
}

void Bspline::setKnotvector(KnotVector<double> &newKnotVector)
{
	knotVector = newKnotVector;
	basisFunctions.setKnotVector(knotVector);
}



int Bspline::findSpanOfValue(double point)
{
	return knotVector.findSpanOfValue(point);
}

std::pair<std::vector<double>, std::vector<double>> Bspline::evaluateAtPoint(double point)
{
	return basisFunctions.evaluateAtPoint(point);
}



void Bspline::plot(int resolution)
{
	double firstKnot = knotVector(0);
	double lastKnot = knotVector(knotVector.getSize()-1);
	std::vector<std::vector<double>> ValuesOfBasisFunctions;
	std::vector<double> evaluationPoints;
	for (int i = 0; i < resolution; i++)
	{
		double currentStep = firstKnot + (double)(i) * ((lastKnot - firstKnot) / ((double)(resolution - 1)));
		evaluationPoints.push_back(currentStep);
		std::vector<double> values = basisFunctions.evaluateAtPoint(currentStep, true).first;
		ValuesOfBasisFunctions.push_back(values);
	}

	std::string filename("basis.dat");
	std::ofstream plotBspline(filename);
	plotBspline << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
	plotBspline << "zone t= " << "\"1\"" << ",i=" << resolution * ValuesOfBasisFunctions[0].size() << ",j=" << resolution << "\n";

	for (int i = 0; i < ValuesOfBasisFunctions[0].size(); i++)
	{
		for (int j = 0; j < resolution; j++)
		{
			plotBspline << evaluationPoints[j] << " " << ValuesOfBasisFunctions[j][i] << "\n";
		}
	}
	plotBspline.close();
}
