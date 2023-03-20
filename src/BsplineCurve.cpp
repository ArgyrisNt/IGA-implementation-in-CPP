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
	plotCurve << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
   	plotCurve << "zone t= " << "\"1\"" << ",i=" << resolution + 1 << ",j=" << resolution + 1 << "\n";

	for (int i = 0; i < steps.size(); ++i)
	{
		int span = multiBspline.findSpanOfValue(steps[i], 0);

		Vertex<double> coordinates = evaluateAtPoint(steps[i]);
		plotCurve << coordinates.x << " " << coordinates.y << "\n";
	}
	plotCurve.close();
}

void BsplineCurve::plotVectorOnEntity(const int resolution, const std::vector<double> &zCoordinate, const std::string &filename)
{
	std::vector<double> steps = multiBspline.getKnotvector(0).linspace(resolution);

	std::ofstream plotCurve(filename);
	plotCurve << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
   	plotCurve << "zone t= " << "\"1\"" << ",i=" << resolution + 1 << ",j=" << resolution + 1 << "\n";

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
		plotCurve << coord_x << " " << coord_z << "\n";
	}
	plotCurve.close();
}
