#ifndef H_TRIMMINGCURVE
#define H_TRIMMINGCURVE

#include <iostream>
#include <vector>

class TrimmingCurve
{
public:
    TrimmingCurve() {}
    TrimmingCurve(std::pair<double, double> &_center, double _radius) : center(_center), radius(_radius) {}
    TrimmingCurve(std::pair<double, double> &&_center, double _radius) : center(_center), radius(_radius) {}

    void plot();
    double projectionOfPoint(std::pair<double, double>& point);
    std::pair<double, double> evaluate(double t);
    std::pair<double, double> evaluateDerivative(double t);
    bool isPointOutside(std::pair<double, double> &point);
    double find_s_given_t(double t, double minimum, double maximum);
    double find_t_given_s(double s, double minimum, double maximum);

    std::pair<double, double> center;
    double radius;
};

#include "..\src\TrimmingCurve.cpp"

#endif