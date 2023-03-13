#include "..\include\BsplineSurface.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

BsplineSurface &BsplineSurface::operator=(const BsplineSurface &old)
{
	bspline_x = old.bspline_x;
	bspline_y = old.bspline_y;
	controlPoints = old.controlPoints;

	return *this;
}



void BsplineSurface::setControlPoints(std::vector<std::vector<double>> &new_controlPoints)
{
    controlPoints = new_controlPoints;
}

void BsplineSurface::setBspline_x(Bspline &new_bspline_x)
{
    bspline_x = new_bspline_x;
}

void BsplineSurface::setBspline_y(Bspline &new_bspline_y)
{
    bspline_y = new_bspline_y;
}

Bspline &BsplineSurface::getBspline_x()
{
    return bspline_x;
}

Bspline &BsplineSurface::getBspline_y()
{
    return bspline_y;
}

std::vector<std::vector<double>> &BsplineSurface::getControlPoints()
{
    return controlPoints;
}



Vertex<double> BsplineSurface::evaluateAtPoint(Vertex<double> &&point)
{
	int span_i = bspline_x.findSpanOfValue(point.x);
	std::vector<double> x_ValuesOfbasisFunctions = bspline_x.evaluateAtPoint(point.x).first;
	int span_j = bspline_y.findSpanOfValue(point.y);
	std::vector<double> y_ValuesOfbasisFunctions = bspline_y.evaluateAtPoint(point.y).first;

	double coordinate_x = 0.0, coordinate_y = 0.0;
	for (int ii = 0; ii < x_ValuesOfbasisFunctions.size(); ii++)
	{
		for (int jj = 0; jj < y_ValuesOfbasisFunctions.size(); jj++)
		{
			int x_index = span_i - bspline_x.getDegree() + ii;
			int y_index = span_j - bspline_y.getDegree() + jj;
			int index = x_index * bspline_y.getNumberOfBasisFunctions() + y_index;
			coordinate_x += x_ValuesOfbasisFunctions[ii] * y_ValuesOfbasisFunctions[jj] * controlPoints[index][0];
			coordinate_y += x_ValuesOfbasisFunctions[ii] * y_ValuesOfbasisFunctions[jj] * controlPoints[index][1];
		}
	}

	return Vertex<double>(coordinate_x, coordinate_y);
}



void BsplineSurface::uniformRefine_x()
{
	std::vector<std::vector<std::vector<double>>> all_new_controlPoints;
	std::vector<std::vector<double>> controlPoints_x;
	KnotVector<double> new_knotvector;
	for (int j = 0; j < bspline_y.getNumberOfBasisFunctions(); j++)
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
	new_knotvector.setWeights(new_weights);
	new_knotvector.computeDistinctKnots();
	bspline_x.setKnotvector(new_knotvector);
	(*this).setBspline_x(bspline_x);
	(*this).setControlPoints(controlPoints);
}

void BsplineSurface::uniformRefine_y()
{
	std::vector<std::vector<std::vector<double>>> all_new_controlPoints;
	std::vector<std::vector<double>> controlPoints_y;
	KnotVector<double> new_knotvector;
	for (int j = 0; j < bspline_x.getNumberOfBasisFunctions(); j++)
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
	new_knotvector.setWeights(new_weights);
	new_knotvector.computeDistinctKnots();
	bspline_y.setKnotvector(new_knotvector);
	(*this).setBspline_y(bspline_y);
	(*this).setControlPoints(controlPoints);
}

std::vector<std::vector<double>> BsplineSurface::pointsOfParametricCurve(int direction, int level)
{
	std::vector<std::vector<double>> controlPointsOnDirection; 
	switch (direction)
	{
	case 0: // x
		for (int i = 0; i < controlPoints.size(); i++)
		{
			if (i == 0 || (i % bspline_y.getNumberOfBasisFunctions()) == 0)
				controlPointsOnDirection.push_back(controlPoints[i + level]);
		}
		break;
	case 1: // y
		for (int i = 0; i < controlPoints.size(); i++)
		{
			if (i >= 0 && (i < bspline_y.getNumberOfBasisFunctions()))
				controlPointsOnDirection.push_back(controlPoints[i + level * bspline_y.getNumberOfBasisFunctions()]);
		}
		break;
	default:
		std::cout << "Invalid dimension. Valid directions are 0,1" << std::endl;
		throw std::invalid_argument("Invalid direction");
		break;
	}

	return controlPointsOnDirection;
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



void BsplineSurface::plot(int resolution)
{
	double firstKnot_x = bspline_x.getKnotvector()(0);
	double lastKnot_x = bspline_x.getKnotvector()(bspline_x.getKnotvector().getSize() - 1);
	double firstKnot_y = bspline_y.getKnotvector()(0);
	double lastKnot_y = bspline_y.getKnotvector()(bspline_y.getKnotvector().getSize() - 1);

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
