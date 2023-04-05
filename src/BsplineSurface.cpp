#include "..\include\BsplineSurface.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

BsplineSurface &BsplineSurface::operator=(const BsplineSurface &old)
{
	controlPoints = old.controlPoints;
	trimmingCurve = old.trimmingCurve;

	return *this;
}


Vertex<double> BsplineSurface::evaluateAtPoint(const Vertex<double> &point)
{
	int span_i = multiBspline.getBspline(0).findSpanOfValue(point.x);
	std::vector<double> x_ValuesOfbasisFunctions = multiBspline.getBspline(0).evaluateAtPoint(point.x).first;
	int span_j = multiBspline.getBspline(1).findSpanOfValue(point.y);
	std::vector<double> y_ValuesOfbasisFunctions = multiBspline.getBspline(1).evaluateAtPoint(point.y).first;

	double coordinate_x = 0.0, coordinate_y = 0.0;
	for (int ii = 0; ii < x_ValuesOfbasisFunctions.size(); ++ii)
	{
		for (int jj = 0; jj < y_ValuesOfbasisFunctions.size(); ++jj)
		{
			int x_index = span_i - multiBspline.getBspline(0).getDegree() + ii;
			int y_index = span_j - multiBspline.getBspline(1).getDegree() + jj;
			int index = x_index * multiBspline.getBspline(1).getNumberOfBasisFunctions() + y_index;
			coordinate_x += x_ValuesOfbasisFunctions[ii] * y_ValuesOfbasisFunctions[jj] * controlPoints[index].x;
			coordinate_y += x_ValuesOfbasisFunctions[ii] * y_ValuesOfbasisFunctions[jj] * controlPoints[index].y;
		}
	}

	return Vertex<double>(coordinate_x, coordinate_y);
}


void BsplineSurface::uniformRefine_x()
{
	std::vector<std::vector<Vertex<double>>> all_new_controlPoints;
	std::vector<Vertex<double>> controlPoints_x;
	KnotVector<double> new_knotvector;
	for (int j = 0; j < multiBspline.getBspline(1).getNumberOfBasisFunctions(); ++j)
	{
		new_knotvector = multiBspline.getBspline(0).getKnotvector();
		controlPoints_x = XparametricCurvePoints(j);
		refineParametricCurve(new_knotvector, controlPoints_x);
		all_new_controlPoints.push_back(controlPoints_x);
	}

	controlPoints = {};
	for (int i = 0; i < all_new_controlPoints[0].size(); ++i)
	{
		for (auto group : all_new_controlPoints)
		{
			controlPoints.push_back(group[i]);
		}
	}

	new_knotvector.computeDistinctKnots();
	multiBspline.getBspline(0).setKnotvector(new_knotvector);
	multiBspline.setBspline(multiBspline.getBspline(0), 0);
	(*this).setControlPoints(controlPoints);
}

void BsplineSurface::uniformRefine_y()
{
	std::vector<std::vector<Vertex<double>>> all_new_controlPoints;
	std::vector<Vertex<double>> controlPoints_y;
	KnotVector<double> new_knotvector;
	for (int j = 0; j < multiBspline.getBspline(0).getNumberOfBasisFunctions(); ++j)
	{
		new_knotvector = multiBspline.getBspline(1).getKnotvector();
		controlPoints_y = YparametricCurvePoints(j);
		refineParametricCurve(new_knotvector, controlPoints_y);
		all_new_controlPoints.push_back(controlPoints_y);
	}

	controlPoints = {};
	for (auto group : all_new_controlPoints)
	{
		for (auto el : group) controlPoints.push_back(el);
	}

	new_knotvector.computeDistinctKnots();
	multiBspline.getBspline(1).setKnotvector(new_knotvector);
	multiBspline.setBspline(multiBspline.getBspline(1), 1);
	(*this).setControlPoints(controlPoints);
}

std::vector<Vertex<double>> BsplineSurface::XparametricCurvePoints(const int level)
{
	std::vector<Vertex<double>> XcontrolPointsOnDirection; 
	for (int i = 0; i < controlPoints.size(); ++i)
	{
		if (i == 0 || (i % multiBspline.getBspline(1).getNumberOfBasisFunctions()) == 0)
			XcontrolPointsOnDirection.push_back(controlPoints[i + level]);
	}
	return XcontrolPointsOnDirection;
}

std::vector<Vertex<double>> BsplineSurface::YparametricCurvePoints(const int level)
{
	std::vector<Vertex<double>> YcontrolPointsOnDirection;
	for (int i = 0; i < controlPoints.size(); ++i)
	{
		if (i >= 0 && (i < multiBspline.getBspline(1).getNumberOfBasisFunctions()))
			YcontrolPointsOnDirection.push_back(controlPoints[i + level * multiBspline.getBspline(1).getNumberOfBasisFunctions()]);
	}
	return YcontrolPointsOnDirection;
}

void BsplineSurface::refineParametricCurve(KnotVector<double> &vector, std::vector<Vertex<double>> &points)
{
	KnotVector<double> new_knotvector = vector;
	for (int i = 0; i < vector.getSize() - 1; ++i)
	{
		double a = vector(i);
		double b = vector(i + 1);
		if (almostEqual(a, b)) continue;
		double new_knot = (a + b) / 2.0;
		knotInsertion(new_knotvector, points, new_knot);
	}
	vector = new_knotvector;
}


void BsplineSurface::plot(const int resolution, const std::string &filename)
{
	std::vector<double> i_steps = multiBspline.getBspline(0).getKnotvector().linspace(resolution);
	std::vector<double> j_steps = multiBspline.getBspline(1).getKnotvector().linspace(resolution);

	std::ofstream plotSurface(filename);
	plotSurface << "variables= " << "\"x\"" << "," << "\"y\""<< "\n";
	plotSurface << "zone t= "<< "\"1\"" << ",i=" << resolution + 1 << ",j=" << resolution + 1<< "\n";

	for (int i = 0; i < i_steps.size(); ++i)
	{
		for (int j = 0; j < j_steps.size(); ++j)
		{
			if (trimmingCurve.isCartesianPointInside(i_steps[i], j_steps[j])) continue;

			Vertex<double> coordinates = evaluateAtPoint(Vertex<double>(i_steps[i], j_steps[j]));

			plotSurface << coordinates.x << " " << coordinates.y << "\n";
		}
	}
	plotSurface.close();
	plotControlPoints("controlPoints.csv");
}

void BsplineSurface::plotVectorOnEntity(const int resolution, const std::vector<double> &zCoordinate, const std::string &filename)
{
	std::vector<double> i_steps = multiBspline.getBspline(0).getKnotvector().linspace(resolution);
	std::vector<double> j_steps = multiBspline.getBspline(1).getKnotvector().linspace(resolution);

	std::ofstream plotSurface(filename);
	plotSurface << "variables= " << "\"x\"" << "," << "\"y\"" << "," << "\"sol\"" << "\n";
    plotSurface << "zone t= " << "\"1\"" << ",i=" << resolution+1 << ",j=" << resolution+1 << "\n";

	for (int i = 0; i < i_steps.size(); ++i)
	{
		for (int j = 0; j < j_steps.size(); ++j)
		{
			if (trimmingCurve.isCartesianPointInside(i_steps[i], j_steps[j])) continue;

			int span_i = multiBspline.getBspline(0).findSpanOfValue(i_steps[i]);
			std::vector<double> bVal_i = multiBspline.getBspline(0).evaluateAtPoint(i_steps[i]).first;
			int span_j = multiBspline.getBspline(1).findSpanOfValue(j_steps[j]);
			std::vector<double> bVal_j = multiBspline.getBspline(1).evaluateAtPoint(j_steps[j]).first;

			double coord_x = 0.0, coord_y = 0.0, coord_z = 0.0;
			for (int kkx = 0; kkx < bVal_i.size(); ++kkx)
			{
				for (int kky = 0; kky < bVal_j.size(); ++kky)
				{
					int i1 = span_i - multiBspline.getBspline(0).getDegree() + kkx;
					int i2 = span_j - multiBspline.getBspline(1).getDegree() + kky;
					int my = i1 * multiBspline.getBspline(1).getNumberOfBasisFunctions() + i2;
					coord_x += bVal_i[kkx] * bVal_j[kky] * controlPoints[my].x;
					coord_y += bVal_i[kkx] * bVal_j[kky] * controlPoints[my].y;
					coord_z += bVal_i[kkx] * bVal_j[kky] * zCoordinate[my];
				}
			}
			plotSurface << coord_x << " " << coord_y << " " << coord_z << "\n";
		}
	}
	plotSurface.close();
}