#include <iostream>
#include <fstream>
#include "math.h"
#include "..\include\TrimmingCurve.h"


void TrimmingCurve::plot()
{
	std::string filename("trimming_curve.dat");
	std::ofstream plotTrimmedCurve(filename);
	plotTrimmedCurve << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
    plotTrimmedCurve << "zone t= " << "\"1\"" << ",i=" << 63 << ",j=" << 63 << "\n";
	double u = 0.0;
	while (u <= 2 * 3.14159265)
	{
        Vertex<double> point = evaluate(u);
        plotTrimmedCurve << point.x << " " << point.y << "\n";
		u += 0.01;
	}
	plotTrimmedCurve.close();
}

Vertex<double> TrimmingCurve::evaluate(double t)
{
    assert(t >= 0.0 && t < 2.0 * 3.14159265);
    return Vertex<double>(center.x + radius * cos(t), center.y + radius * sin(t));
}

Vertex<double> TrimmingCurve::evaluateDerivative(double t)
{
    assert(t >= 0.0 && t < 2.0 * 3.14159265);
    return Vertex<double>(-radius * sin(t), radius * cos(t));
}

double TrimmingCurve::projectionOfPoint(Vertex<double> &point)
{
    std::vector<double> possible_projs;
    double u = 0.0;
    while (u <= 2 * 3.14159265)
    {
        Vertex<double> tangent_vec = evaluateDerivative(u);
        Vertex<double> di(point.x - evaluate(u).x, point.y - evaluate(u).y);
        double dot = tangent_vec.x * di.x + tangent_vec.y * di.y;
        double d = sqrt(di.x * di.x + di.y * di.y);
        if (dot <= 1e-3 && dot >= -1e-3)
        {
            possible_projs.push_back(d);
        }
        u += 0.001;
    }

    double min = possible_projs[0];
    for (int i = 1; i < possible_projs.size(); i++)
    {
        if (possible_projs[i] < min)
            min = possible_projs[i];
    }

    return min;
}

bool TrimmingCurve::isPointOutside(Vertex<double> &point)
{
    std::vector<double> possible_projs;
    std::vector<Vertex<double>> possible_tangents;
    std::vector<Vertex<double>> possible_dis;
    double u = 0.0;
    while (u <= 2 * 3.14159265)
    {
        Vertex<double> tangent_vec = evaluateDerivative(u);
        Vertex<double> di(point.x - evaluate(u).x, point.y - evaluate(u).y);
        double dot = tangent_vec.x * di.x + tangent_vec.y * di.y;
        double d = sqrt(di.x * di.x + di.y * di.y);
        if (dot <= 1e-3 && dot >= -1e-3)
        {
            possible_projs.push_back(d);
            possible_tangents.push_back(tangent_vec);
            possible_dis.push_back(di);
        }
        u += 0.001;
    }

    double min = 0;
    for (int i = 1; i < possible_projs.size(); i++)
    {
        if (possible_projs[i] < possible_projs[min])
            min = i;
    }

    std::vector<double> v1{possible_dis[min].x, possible_dis[min].y, 0.0};
    std::vector<double> v2{possible_tangents[min].x, possible_tangents[min].y, 0.0};
    double w1 = v1[1] * v2[2] - v1[2] * v2[1];
    double w2 = v1[2] * v2[0] - v1[0] * v2[2];
    double w3 = v1[0] * v2[1] - v1[1] * v2[0];

    if (w3 < 0 && w3 < -1e-5) return true;
    else return false;
}

double TrimmingCurve::find_s_given_t(double t, double minimum, double maximum)
{
    double s = 0.0;
    double u = 0.0;
    while (u <= 2 * 3.14159265)
    {
        Vertex<double> point = evaluate(u);
        if (fabs(point.y - t) < 1e-3)
        {
            if (point.x >= minimum && point.x <= maximum)
            {
                s = point.x;
                break;
            }
        }
        u += 0.001;
    }
    return s;
}

double TrimmingCurve::find_t_given_s(double s, double minimum, double maximum)
{
    double t = 0.0;
    double u = 0.0;
    while (u <= 2 * 3.14159265)
    {
        Vertex<double> point = evaluate(u);
        if (fabs(point.x - s) < 1e-3)
        {
            if (point.y >= minimum && point.y <= maximum)
            {
                t = point.y;
                break;
            }
        }
        u += 0.001;
    }
    return t;
}

bool TrimmingCurve::isCartesianPointInside(double x, double y)
{
    return (std::pow(x - center.x, 2) + std::pow(y - center.y, 2)) < std::pow(radius, 2);
}