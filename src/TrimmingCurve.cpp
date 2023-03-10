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
        std::pair<double, double> point = evaluate(u);
        plotTrimmedCurve << point.first << " " << point.second << "\n";
		u += 0.01;
	}
	plotTrimmedCurve.close();
}

std::pair<double, double> TrimmingCurve::evaluate(double t)
{
    if (t < 0.0 || t >= 2.0 * 3.14159265)
        throw std::invalid_argument("Invalid parameter value");
    return std::make_pair(center.first + radius * cos(t), center.second + radius * sin(t));
}

std::pair<double, double> TrimmingCurve::evaluateDerivative(double t)
{
    if (t < 0.0 || t >= 2.0 * 3.14159265)
        throw std::invalid_argument("Invalid parameter value");
    return std::make_pair(-radius * sin(t), radius * cos(t));
}

double TrimmingCurve::projectionOfPoint(std::pair<double, double>& point)
{
    std::vector<double> possible_projs;
    double u = 0.0;
    while (u <= 2 * 3.14159265)
    {
        std::pair<double, double> tangent_vec = evaluateDerivative(u);
        std::pair<double, double> di = std::make_pair(point.first - evaluate(u).first, point.second - evaluate(u).second);
        double dot = tangent_vec.first * di.first + tangent_vec.second * di.second;
        double d = sqrt(di.first * di.first + di.second * di.second);
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

bool TrimmingCurve::isPointOutside(std::pair<double, double>& point)
{
    std::vector<double> possible_projs;
    std::vector<std::pair<double, double>> possible_tangents;
    std::vector<std::pair<double, double>> possible_dis;
    double u = 0.0;
    while (u <= 2 * 3.14159265)
    {
        std::pair<double, double> tangent_vec = evaluateDerivative(u);
        std::pair<double, double> di = std::make_pair(point.first - evaluate(u).first, point.second - evaluate(u).second);
        double dot = tangent_vec.first * di.first + tangent_vec.second * di.second;
        double d = sqrt(di.first * di.first + di.second * di.second);
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

    std::vector<double> v1{possible_dis[min].first, possible_dis[min].second, 0.0};
    std::vector<double> v2{possible_tangents[min].first, possible_tangents[min].second, 0.0};
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
        std::pair<double, double> point = evaluate(u);
        if (fabs(point.second - t) < 1e-3)
        {
            if (point.first >= minimum && point.first <= maximum)
            {
                s = point.first;
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
        std::pair<double, double> point = evaluate(u);
        if (fabs(point.first - s) < 1e-3)
        {
            if (point.second >= minimum && point.second <= maximum)
            {
                t = point.second;
                break;
            }
        }
        u += 0.001;
    }
    return t;
}