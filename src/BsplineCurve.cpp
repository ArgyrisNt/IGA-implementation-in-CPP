#include "..\include\BsplineCurve.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

Bspline &BsplineCurve::getBspline_x()
{
    return bspline_x;
}

std::vector<std::vector<double>> &BsplineCurve::getControlPoints()
{
    return controlPoints;
}



std::pair<double, double> BsplineCurve::evaluateAtPoint(double point)
{
	int span = bspline_x.findSpanOfValue(point);
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
