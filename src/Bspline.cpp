#include "..\include\Bspline.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

Bspline::Bspline(int newDegree, KnotVector<double>& newKnotVector, std::vector<double> &newWeights)
{
	weights = newWeights;
	knotVector = newKnotVector;
	degree = newDegree;
	numberOfBasisFunctions = knotVector.getSize() - degree - 1;
}

Bspline::Bspline(const Bspline &bspline)
{
	knotVector = bspline.knotVector;
	weights = bspline.weights;
	degree = bspline.degree;
	numberOfBasisFunctions = bspline.numberOfBasisFunctions;
}

Bspline::~Bspline() {}

Bspline &Bspline::operator=(const Bspline &bspline)
{
	knotVector = bspline.knotVector;
	weights = bspline.weights;
	degree = bspline.degree;
	numberOfBasisFunctions = bspline.numberOfBasisFunctions;

	return *this;
}

void Bspline::setWeights(std::vector<double> new_weights)
{
	for (auto weight : new_weights)
	{
		if (weight < 0 || weight > 1)
		{
			std::cout << "Invalid weights. Valid weights are within [0,1]." << std::endl;
			throw std::invalid_argument("Invalid weights");
			break;
		}
	}

	weights = new_weights;
}

void Bspline::basisFunctionsOfDegree(int level, double value, std::vector<double> &valuesOfBasisFunctions, std::vector<double>& derivativesOfBasisFunctions)
{
	std::vector<double> oldValuesOfBasisFunctions = valuesOfBasisFunctions;
	for (int j = 0; j < numberOfBasisFunctions; j++)
	{
		double firstCoefficient = 0.0, secondCoefficient = 0.0;

		double denominator = knotVector(j + level) - knotVector(j);
		if (denominator != 0.0) firstCoefficient = (value - knotVector(j)) / denominator;

		denominator = knotVector(j + level + 1) - knotVector(j + 1);
		if (denominator != 0.0) secondCoefficient = (knotVector(j + level + 1) - value) / denominator;

		if (j == knotVector.getSize() - degree - 2) valuesOfBasisFunctions[j] = firstCoefficient * oldValuesOfBasisFunctions[j];
		else valuesOfBasisFunctions[j] = firstCoefficient * oldValuesOfBasisFunctions[j] + secondCoefficient * oldValuesOfBasisFunctions[j + 1];

		if (level == degree)
		{
			denominator = knotVector(j + level) - knotVector(j);
			if (denominator != 0.0) firstCoefficient = oldValuesOfBasisFunctions[j] / (knotVector(j + level) - knotVector(j));
			denominator = knotVector(j + level + 1) - knotVector(j + 1);
			if (denominator != 0.0) secondCoefficient = oldValuesOfBasisFunctions[j + 1] / (knotVector(j + level + 1) - knotVector(j + 1));
			derivativesOfBasisFunctions[j] = degree * (firstCoefficient - secondCoefficient);
		}
	}
}

void Bspline::computeActiveBasisFunctions(double value, std::vector<double> &valuesOfBasisFunctions, std::vector<double> &derivativesOfBasisFunctions)
{
	// if value in [ui,ui+1), then active functions are Ni-p,...,Ni
	std::vector<double> activeValuesOfBasisFunctions;
	std::vector<double> activeDerivativesOfBasisFunctions;
	int span = knotVector.findSpanOfValue(value);
	double sumValues = 0.0, sumDerivatives = 0.0;
	for (int i = span - degree; i <= span; i++)
	{
		activeValuesOfBasisFunctions.push_back(valuesOfBasisFunctions[i]);
		activeDerivativesOfBasisFunctions.push_back(derivativesOfBasisFunctions[i]);
		sumValues += valuesOfBasisFunctions[i] * weights[i];
		sumDerivatives += derivativesOfBasisFunctions[i] * weights[i];
	}

	std::vector<double> activeNurbs, activeNurbsDerivatives;
	for (int i = 0; i < activeValuesOfBasisFunctions.size(); i++)
	{
		double weightedBasis = activeValuesOfBasisFunctions[i] * weights[i];
		double weightedDerivatives = activeDerivativesOfBasisFunctions[i] * weights[i];
		activeNurbs.push_back(weightedBasis / sumValues);
		activeNurbsDerivatives.push_back((weightedDerivatives * sumValues - weightedBasis * sumDerivatives) / (sumValues * sumValues));
	}

	valuesOfBasisFunctions = activeNurbs;
	derivativesOfBasisFunctions = activeNurbsDerivatives;
}

std::vector<double> Bspline::initializeBasisFunctions(double value)
{
	std::vector<double> valuesOfBasisFunctions(numberOfBasisFunctions, 0.0);
	for (int j = 0; j < numberOfBasisFunctions; j++)
	{
		bool isBetweenTheseTwoKnots = (value >= knotVector(j) && value < knotVector(j + 1));
		if (isBetweenTheseTwoKnots) valuesOfBasisFunctions[j] = 1.0;
	}

	bool isLastKnot = (value == knotVector(knotVector.getSize() - 1));
	if (isLastKnot) valuesOfBasisFunctions[numberOfBasisFunctions - 1] = 1.0;

	return valuesOfBasisFunctions;
}

std::pair<std::vector<double>, std::vector<double>> Bspline::evaluateAtPoint(const double value, bool allBasisFunctions)
{
	std::vector<double> valuesOfBasisFunctions = initializeBasisFunctions(value);
	std::vector<double> derivativesOfBasisFunctions(numberOfBasisFunctions);

	for (int i = 1; i < degree + 1; i++)
	{
		basisFunctionsOfDegree(i, value, valuesOfBasisFunctions, derivativesOfBasisFunctions);		
	}

	if (!allBasisFunctions)
	{
		computeActiveBasisFunctions(value, valuesOfBasisFunctions, derivativesOfBasisFunctions);
	}

	return std::make_pair(valuesOfBasisFunctions, derivativesOfBasisFunctions);
}

void Bspline::plot(int resolution)
{
	double firstKnot = getKnotvector()(0);
	double lastKnot = getKnotvector()(getKnotvector().getSize() - 1);
	std::vector<std::vector<double>> ValuesOfBasisFunctions;
	std::vector<double> evaluationPoints;
	for (int i = 0; i < resolution; i++)
	{
		double currentStep = firstKnot + (double)(i) * ((lastKnot - firstKnot) / ((double)(resolution - 1)));
		evaluationPoints.push_back(currentStep);
		std::vector<double> values = evaluateAtPoint(currentStep, true).first;
		ValuesOfBasisFunctions.push_back(values);
	}

	std::string filename("basis.dat");
	std::ofstream plotBspline(filename);
	plotBspline << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
	plotBspline << "zone t= " << "\"1\"" << ",i=" << ValuesOfBasisFunctions.size() * ValuesOfBasisFunctions[0].size() << ",j=" << ValuesOfBasisFunctions.size() << "\n";

	for (int i = 0; i < ValuesOfBasisFunctions[0].size(); i++)
	{
		for (int j = 0; j < ValuesOfBasisFunctions.size(); j++)
		{
			plotBspline << evaluationPoints[j] << " " << ValuesOfBasisFunctions[j][i] << "\n";
		}
	}
	plotBspline.close();
}

std::pair<double, double> BsplineCurve::evaluateAtPoint(double point)
{
	int span = bspline_x.getKnotvector().findSpanOfValue(point);
	std::vector<double> ValuesOfbasisFunctions = bspline_x.evaluateAtPoint(point).first;

	double coordinate_x = 0.0, coordinate_y = 0.0;
	for (int j = 0; j < ValuesOfbasisFunctions.size(); j++)
	{
		coordinate_x += ValuesOfbasisFunctions[j] * controlPoints[span - bspline_x.getDegree() + j][0];
		coordinate_y += ValuesOfbasisFunctions[j] * controlPoints[span - bspline_x.getDegree() + j][1];
	}

	return std::make_pair(coordinate_x, coordinate_y);
}

void BsplineCurve::plot(int resolution)
{
	double firstKnot = bspline_x.getKnotvector()(0);
	double lastKnot = bspline_x.getKnotvector()(bspline_x.getKnotvector().getSize() - 1);

	std::string filename("curve.dat");
	std::ofstream plotCurve(filename);
	plotCurve << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
   	plotCurve << "zone t= " << "\"1\"" << ",i=" << resolution << ",j=" << resolution << "\n";

	for (int i = (int)(firstKnot); i < resolution; i++)
	{
		double currentStep = firstKnot + (double)(i) * ((lastKnot - firstKnot) / ((double) (resolution - 1)));
		std::pair<double,double> coordinates = evaluateAtPoint(currentStep);
		plotCurve << coordinates.first << " " << coordinates.second << "\n";
	}
	plotCurve.close();
}

BsplineSurface &BsplineSurface::operator=(const BsplineSurface &old)
{
	bspline_x = old.bspline_x;
	bspline_y = old.bspline_y;
	controlPoints = old.controlPoints;

	return *this;
}

Vertex<double> BsplineSurface::evaluateAtPoint(Vertex<double> &&point)
{
	int span_i = getBspline_x().getKnotvector().findSpanOfValue(point.x);
	std::vector<double> x_ValuesOfbasisFunctions = getBspline_x().evaluateAtPoint(point.x).first;
	int span_j = getBspline_y().getKnotvector().findSpanOfValue(point.y);
	std::vector<double> y_ValuesOfbasisFunctions = getBspline_y().evaluateAtPoint(point.y).first;

	double coordinate_x = 0.0, coordinate_y = 0.0;
	for (int ii = 0; ii < x_ValuesOfbasisFunctions.size(); ii++)
	{
		for (int jj = 0; jj < y_ValuesOfbasisFunctions.size(); jj++)
		{
			int x_index = span_i - getBspline_x().getDegree() + ii;
			int y_index = span_j - getBspline_y().getDegree() + jj;
			int index = x_index * getBspline_y().getNumberOfBasisFunctions() + y_index;
			coordinate_x += x_ValuesOfbasisFunctions[ii] * y_ValuesOfbasisFunctions[jj] * getControlPoints()[index][0];
			coordinate_y += x_ValuesOfbasisFunctions[ii] * y_ValuesOfbasisFunctions[jj] * getControlPoints()[index][1];
		}
	}

	return Vertex<double>(coordinate_x, coordinate_y);
}

void BsplineSurface::plot(int resolution)
{
	double firstKnot_x = getBspline_x().getKnotvector()(0);
	double lastKnot_x = getBspline_x().getKnotvector()(getBspline_x().getKnotvector().getSize() - 1);
	double firstKnot_y = getBspline_y().getKnotvector()(0);
	double lastKnot_y = getBspline_y().getKnotvector()(getBspline_y().getKnotvector().getSize() - 1);

	std::string filename("surface.dat");
	std::ofstream plotSurface(filename);
	plotSurface << "variables= " << "\"x\"" << "," << "\"y\""<< "\n";
	plotSurface << "zone t= "<< "\"1\"" << ",i=" << resolution << ",j=" << resolution << "\n";

	for (int j = (int)(firstKnot_y); j < resolution; j++)
	{
		for (int i = (int)(firstKnot_x); i < resolution; i++)
		{
			double currentStep_i = firstKnot_x + (double)(i) * ((lastKnot_x - firstKnot_x) / ((double)(resolution)));
			double currentStep_j = firstKnot_y + (double)(j) * ((lastKnot_y - firstKnot_y) / ((double)(resolution)));
			Vertex<double> coordinates = evaluateAtPoint(Vertex<double>(currentStep_i, currentStep_j));
			plotSurface << coordinates.x << " " << coordinates.y << "\n";
		}
	}
}

void BsplineSurface::knotInsertion(KnotVector<double>& vector, std::vector<std::vector<double>> &points, double newKnot)
{
	KnotVector<double> new_knotvector = vector;
	int span = vector.findSpanOfValue(newKnot);
	new_knotvector.insert(span + 1, newKnot);

	std::vector<std::vector<double>> newPoints;
	for (int j = span - vector.getDegree() + 1; j <= span; j++)
	{
		double alpha = (newKnot - vector(j)) / (vector(j + vector.getDegree()) - vector(j));
		double coord_x = (1.0 - alpha) * points[j - 1][0] + alpha * points[j][0];
		double coord_y = (1.0 - alpha) * points[j - 1][1] + alpha * points[j][1];
		newPoints.push_back({coord_x, coord_y});
	}

	std::vector<std::vector<double>> FinalNewPoints;
	for (int i = 0; i <= span - vector.getDegree(); i++)
	{
		FinalNewPoints.push_back(points[i]);
	}
	for (int ii = 0; ii < newPoints.size(); ii++)
	{
		FinalNewPoints.push_back(newPoints[ii]);
	}
	for (int ii = span; ii < points.size(); ii++)
	{
		FinalNewPoints.push_back(points[ii]);
	}

	points = FinalNewPoints;
	vector = new_knotvector;
}

void BsplineSurface::refineParametricCurve(KnotVector<double> &vector, std::vector<std::vector<double>> &points)
{
	KnotVector<double> new_knotvector = vector;
	for (int i = 0; i < vector.getSize() - 1; i++)
	{
		double a = vector(i);
		double b = vector(i + 1);
		if (a == b) continue;
		double new_knot = (a + b) / 2.0;
		knotInsertion(new_knotvector, points, new_knot);
	}
	vector = new_knotvector;
}

std::vector<std::vector<double>> BsplineSurface::pointsOfParametricCurve(int direction, int level)
{
	std::vector<std::vector<double>> controlPointsOnDirection; 
	switch (direction)
	{
	case 0: // x
		for (int i = 0; i < controlPoints.size(); i++)
		{
			if (i == 0 || (i % bspline_y.getWeights().size()) == 0)
				controlPointsOnDirection.push_back(controlPoints[i + level]);
		}
		break;
	case 1: // y
		for (int i = 0; i < controlPoints.size(); i++)
		{
			if (i >= 0 && (i < bspline_y.getWeights().size()))
				controlPointsOnDirection.push_back(controlPoints[i + level * bspline_y.getWeights().size()]);
		}
		break;
	default:
		std::cout << "Invalid dimension. Valid directions are 0,1" << std::endl;
		throw std::invalid_argument("Invalid direction");
		break;
	}

	return controlPointsOnDirection;
}

void BsplineSurface::uniformRefine_x()
{
	std::vector<std::vector<std::vector<double>>> all_new_controlPoints;
	std::vector<std::vector<double>> controlPoints_x;
	KnotVector<double> new_knotvector;
	for (int j = 0; j < bspline_y.getWeights().size(); j++)
	{
		new_knotvector = bspline_x.getKnotvector();
		controlPoints_x = pointsOfParametricCurve(0, j);
		refineParametricCurve(new_knotvector, controlPoints_x);
		all_new_controlPoints.push_back(controlPoints_x);
	}

	controlPoints = {};
	for (int i = 0; i < all_new_controlPoints[0].size(); i++)
	{
		for (auto group : all_new_controlPoints)
		{
			controlPoints.push_back(group[i]);
		}
	}

	std::vector<double> new_weights(new_knotvector.getSize() - bspline_x.getDegree() - 1, 1.0);
	Bspline new_bspline(bspline_x.getDegree(), new_knotvector, new_weights);
	(*this).setBspline_x(new_bspline);
	(*this).setControlPoints(controlPoints);
}

void BsplineSurface::uniformRefine_y()
{
	std::vector<std::vector<std::vector<double>>> all_new_controlPoints;
	std::vector<std::vector<double>> controlPoints_y;
	KnotVector<double> new_knotvector;
	for (int j = 0; j < bspline_x.getWeights().size(); j++)
	{
		new_knotvector = bspline_y.getKnotvector();
		controlPoints_y = pointsOfParametricCurve(1, j);
		refineParametricCurve(new_knotvector, controlPoints_y);
		all_new_controlPoints.push_back(controlPoints_y);
	}

	controlPoints = {};
	for (auto group : all_new_controlPoints)
	{
		for (auto el : group) controlPoints.push_back(el);
	}

	std::vector<double> new_weights(new_knotvector.getSize() - bspline_y.getDegree() - 1, 1.0);
	Bspline new_bspline(bspline_y.getDegree(), new_knotvector, new_weights);
	(*this).setBspline_y(new_bspline);
	(*this).setControlPoints(controlPoints);
}
