#ifndef H_TRIMMINGCURVE
#define H_TRIMMINGCURVE

#include <iostream>
#include <vector>
#include "..\include\Utilities.h"

class TrimmingCurve
{
public:
    TrimmingCurve() {}
    TrimmingCurve(Vertex<double> &_center, double _radius) : center(_center), radius(_radius) {}
    TrimmingCurve(Vertex<double> &&_center, double _radius) : center(_center), radius(_radius) {}

    ~TrimmingCurve() {}

    Vertex<double> evaluate(double t);
    Vertex<double> evaluateDerivative(double t);

    double find_s_given_t(double t, double minimum, double maximum);
    double find_t_given_s(double s, double minimum, double maximum);

    double projectionOfPoint(Vertex<double>& point);
    bool isPointOutside(Vertex<double> &point);
    bool isCartesianPointInside(double x, double y);

    void plot();

    Vertex<double> center;
    double radius;
};

#include "..\src\TrimmingCurve.cpp"

#endif