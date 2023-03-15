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



void BsplineCurve::plot2D(int resolution, std::string filename)
{
	std::vector<double> steps = bspline_x.getKnotvector().linspace(resolution);

	std::ofstream plotCurve(filename);
	plotCurve << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
   	plotCurve << "zone t= " << "\"1\"" << ",i=" << resolution + 1 << ",j=" << resolution + 1 << "\n";

	for (int i = 0; i < steps.size(); i++)
	{
		int span = bspline_x.findSpanOfValue(steps[i]);
		std::vector<double> bVal = bspline_x.evaluateAtPoint(steps[i]).first;

		std::pair<double, double> coordinates = evaluateAtPoint(steps[i]);
		plotCurve << coordinates.first << " " << coordinates.second << "\n";
	}
	plotCurve.close();
}

void BsplineCurve::plot3D(int resolution, std::vector<double> &zCoordinate, std::string filename)
{
	std::vector<double> steps = bspline_x.getKnotvector().linspace(resolution);

	std::ofstream plotCurve(filename);
	plotCurve << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
   	plotCurve << "zone t= " << "\"1\"" << ",i=" << resolution + 1 << ",j=" << resolution + 1 << "\n";

	for (int i = 0; i < steps.size(); i++)
	{
		int span = bspline_x.findSpanOfValue(steps[i]);
		std::vector<double> bVal = bspline_x.evaluateAtPoint(steps[i]).first;

		double coord_x = 0.0, coord_z = 0.0;
		for (int kk = 0; kk < bVal.size(); kk++)
		{
			coord_x += bVal[kk] * controlPoints[span - bspline_x.getDegree() + kk][0];
			coord_z += bVal[kk] * zCoordinate[span - bspline_x.getDegree() + kk];
		}
		plotCurve << coord_x << " " << coord_z << "\n";
	}
	plotCurve.close();
}
