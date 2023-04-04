#include "..\include\BsplineCurve.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>


Vertex<double> BsplineCurve::evaluateAtPoint(const double point)
{
	int span = multiBspline.findSpanOfValue(point, 0);
	std::vector<double> ValuesOfbasisFunctions = multiBspline.getBspline(0).evaluateAtPoint(point).first;

	double coordinate_x = 0.0, coordinate_y = 0.0;
	for (int j = 0; j < ValuesOfbasisFunctions.size(); ++j)
	{
		coordinate_x += ValuesOfbasisFunctions[j] * controlPoints[span - multiBspline.getDegree(0) + j].x;
		coordinate_y += ValuesOfbasisFunctions[j] * controlPoints[span - multiBspline.getDegree(0) + j].y;
	}

	return Vertex<double>(coordinate_x, coordinate_y);
}



void BsplineCurve::plot(const int resolution, const std::string &filename)
{
	std::vector<double> steps = multiBspline.getKnotvector(0).linspace(resolution);

	std::ofstream plotCurve(filename);
	plotCurve << "X,Y\n";
	for (int i = 0; i < steps.size(); ++i)
	{
		int span = multiBspline.findSpanOfValue(steps[i], 0);

		Vertex<double> coordinates = evaluateAtPoint(steps[i]);
		plotCurve << coordinates.x << "," << coordinates.y << "\n";
	}
	plotCurve.close();
	plotControlPoints("controlPoints.csv");
}

void BsplineCurve::plotControlPoints(const std::string &filename)
{
	std::ofstream plotCtrlPts(filename);
	plotCtrlPts << "X,Y\n";
	for (auto it = controlPoints.begin(); it != controlPoints.end(); ++it)
	{
		plotCtrlPts << (*it).x << "," << (*it).y <<"\n";
	}
	plotCtrlPts.close();
}

void BsplineCurve::plotVectorOnEntity(const int resolution, const std::vector<double> &zCoordinate, const std::string &filename)
{
	std::vector<double> steps = multiBspline.getKnotvector(0).linspace(resolution);

	std::ofstream plotCurve(filename);
	plotCurve << "X,Y\n";

	for (int i = 0; i < steps.size(); ++i)
	{
		int span = multiBspline.findSpanOfValue(steps[i], 0);
		std::vector<double> bVal = multiBspline.getBspline(0).evaluateAtPoint(steps[i]).first;

		double coord_x = 0.0, coord_z = 0.0;
		for (int kk = 0; kk < bVal.size(); ++kk)
		{
			coord_x += bVal[kk] * controlPoints[span - multiBspline.getDegree(0) + kk].x;
			coord_z += bVal[kk] * zCoordinate[span - multiBspline.getDegree(0) + kk];
		}
		plotCurve << coord_x << "," << coord_z << "\n";
	}
	plotCurve.close();
}



void BsplineCurve::uniformRefine_x()
{
	KnotVector<double> vector = multiBspline.getBspline(0).getKnotvector();
	KnotVector<double> new_knotvector = vector;
	for (int i = 0; i < vector.getSize() - 1; ++i)
	{
		double a = vector(i);
		double b = vector(i + 1);
		if (almostEqual(a, b)) continue;
		double new_knot = (a + b) / 2.0;
		knotInsertion(new_knotvector, controlPoints, new_knot);
	}
	vector = new_knotvector;

	std::vector<double> new_weights(new_knotvector.getSize() - multiBspline.getBspline(0).getDegree() - 1, 1.0);
	new_knotvector.setWeights(new_weights);
	new_knotvector.computeDistinctKnots();
	multiBspline.getBspline(0).setKnotvector(new_knotvector);
	multiBspline.setBspline(multiBspline.getBspline(0), 0);
	(*this).setControlPoints(controlPoints);
}