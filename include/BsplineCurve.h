#ifndef H_BSPLINECURVE
#define H_BSPLINECURVE

#include "..\include\Bspline.h"
#include <iostream>
#include <vector>

class BsplineCurve
{
public:
    BsplineCurve(Bspline &new_bspline_x, std::vector<std::vector<double>> &new_controlPoints)
        : bspline_x(new_bspline_x), controlPoints(new_controlPoints) {}

    ~BsplineCurve() {}

    std::pair<double, double> evaluateAtPoint(double point);
    void plot2D(int resolution, std::string filename);
    void plot3D(int resolution, std::vector<double> &zCoordinate, std::string filename);

    Bspline &getBspline_x();
    std::vector<std::vector<double>> &getControlPoints();

private:
    Bspline bspline_x;
    std::vector<std::vector<double>> controlPoints;
};

#include "..\src\BsplineCurve.cpp"

#endif