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



void BsplineSurface::setControlPoints(std::vector<std::vector<double>> &new_controlPoints)
{
    controlPoints = new_controlPoints;
}


std::vector<std::vector<double>> &BsplineSurface::getControlPoints()
{
    return controlPoints;
}



Vertex<double> BsplineSurface::evaluateAtPoint(Vertex<double> &&point)
{
	int span_i = getBspline_x().findSpanOfValue(point.x);
	std::vector<double> x_ValuesOfbasisFunctions = getBspline_x().evaluateAtPoint(point.x).first;
	int span_j = getBspline_y().findSpanOfValue(point.y);
	std::vector<double> y_ValuesOfbasisFunctions = getBspline_y().evaluateAtPoint(point.y).first;

	double coordinate_x = 0.0, coordinate_y = 0.0;
	for (int ii = 0; ii < x_ValuesOfbasisFunctions.size(); ii++)
	{
		for (int jj = 0; jj < y_ValuesOfbasisFunctions.size(); jj++)
		{
			int x_index = span_i - getBspline_x().getDegree() + ii;
			int y_index = span_j - getBspline_y().getDegree() + jj;
			int index = x_index * getBspline_y().getNumberOfBasisFunctions() + y_index;
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
	for (int j = 0; j < getBspline_y().getNumberOfBasisFunctions(); j++)
	{
		new_knotvector = getBspline_x().getKnotvector();
		controlPoints_x = XparametricCurvePoints(j);
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

	std::vector<double> new_weights(new_knotvector.getSize() - getBspline_x().getDegree() - 1, 1.0);
	new_knotvector.setWeights(new_weights);
	new_knotvector.computeDistinctKnots();
	getBspline_x().setKnotvector(new_knotvector);
	(*this).setBspline_x(getBspline_x());
	(*this).setControlPoints(controlPoints);
}

void BsplineSurface::uniformRefine_y()
{
	std::vector<std::vector<std::vector<double>>> all_new_controlPoints;
	std::vector<std::vector<double>> controlPoints_y;
	KnotVector<double> new_knotvector;
	for (int j = 0; j < getBspline_x().getNumberOfBasisFunctions(); j++)
	{
		new_knotvector = getBspline_y().getKnotvector();
		controlPoints_y = YparametricCurvePoints(j);
		refineParametricCurve(new_knotvector, controlPoints_y);
		all_new_controlPoints.push_back(controlPoints_y);
	}

	controlPoints = {};
	for (auto group : all_new_controlPoints)
	{
		for (auto el : group) controlPoints.push_back(el);
	}

	std::vector<double> new_weights(new_knotvector.getSize() - getBspline_y().getDegree() - 1, 1.0);
	new_knotvector.setWeights(new_weights);
	new_knotvector.computeDistinctKnots();
	getBspline_y().setKnotvector(new_knotvector);
	(*this).setBspline_y(getBspline_y());
	(*this).setControlPoints(controlPoints);
}

std::vector<std::vector<double>> BsplineSurface::XparametricCurvePoints(int level)
{
	std::vector<std::vector<double>> XcontrolPointsOnDirection; 
	for (int i = 0; i < controlPoints.size(); i++)
	{
		if (i == 0 || (i % getBspline_y().getNumberOfBasisFunctions()) == 0)
			XcontrolPointsOnDirection.push_back(controlPoints[i + level]);
	}
	return XcontrolPointsOnDirection;
}

std::vector<std::vector<double>> BsplineSurface::YparametricCurvePoints(int level)
{
	std::vector<std::vector<double>> YcontrolPointsOnDirection;
	for (int i = 0; i < controlPoints.size(); i++)
	{
		if (i >= 0 && (i < getBspline_y().getNumberOfBasisFunctions()))
			YcontrolPointsOnDirection.push_back(controlPoints[i + level * getBspline_y().getNumberOfBasisFunctions()]);
	}
	return YcontrolPointsOnDirection;
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



void BsplineSurface::plot2D(int resolution, std::string filename)
{
	std::vector<double> i_steps = getBspline_x().getKnotvector().linspace(resolution);
	std::vector<double> j_steps = getBspline_y().getKnotvector().linspace(resolution);

	std::ofstream plotSurface(filename);
	plotSurface << "variables= " << "\"x\"" << "," << "\"y\""<< "\n";
	plotSurface << "zone t= "<< "\"1\"" << ",i=" << resolution + 1 << ",j=" << resolution + 1<< "\n";

	for (int i = 0; i < i_steps.size(); i++)
	{
		for (int j = 0; j < j_steps.size(); j++)
		{
			if (trimmingCurve.isCartesianPointInside(i_steps[i], j_steps[j])) continue;

			Vertex<double> coordinates = evaluateAtPoint(Vertex<double>(i_steps[i], j_steps[j]));

			plotSurface << coordinates.x << " " << coordinates.y << "\n";
		}
	}
	plotSurface.close();
}

void BsplineSurface::plot3D(int resolution, std::vector<double> &zCoordinate, std::string filename)
{
	std::vector<double> i_steps = getBspline_x().getKnotvector().linspace(resolution);
	std::vector<double> j_steps = getBspline_y().getKnotvector().linspace(resolution);

	std::ofstream plotSurface(filename);
	plotSurface << "variables= " << "\"x\"" << "," << "\"y\"" << "," << "\"sol\"" << "\n";
    plotSurface << "zone t= " << "\"1\"" << ",i=" << resolution+1 << ",j=" << resolution+1 << "\n";

	for (int i = 0; i < i_steps.size(); i++)
	{
		for (int j = 0; j < j_steps.size(); j++)
		{
			if (trimmingCurve.isCartesianPointInside(i_steps[i], j_steps[j])) continue;

			int span_i = getBspline_x().findSpanOfValue(i_steps[i]);
			std::vector<double> bVal_i = getBspline_x().evaluateAtPoint(i_steps[i]).first;
			int span_j = getBspline_y().findSpanOfValue(j_steps[j]);
			std::vector<double> bVal_j = getBspline_y().evaluateAtPoint(j_steps[j]).first;

			double coord_x = 0.0, coord_y = 0.0, coord_z = 0.0;
			for (int kkx = 0; kkx < bVal_i.size(); kkx++)
			{
				for (int kky = 0; kky < bVal_j.size(); kky++)
				{
					int i1 = span_i - getBspline_x().getDegree() + kkx;
					int i2 = span_j - getBspline_y().getDegree() + kky;
					int my = i1 * getBspline_y().getNumberOfBasisFunctions() + i2;
					coord_x += bVal_i[kkx] * bVal_j[kky] * controlPoints[my][0];
					coord_y += bVal_i[kkx] * bVal_j[kky] * controlPoints[my][1];
					coord_z += bVal_i[kkx] * bVal_j[kky] * zCoordinate[my];
				}
			}
			plotSurface << coord_x << " " << coord_y << " " << coord_z << "\n";
		}
	}
	plotSurface.close();
}