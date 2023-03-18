#ifndef H_BSPLINECURVE
#define H_BSPLINECURVE

#include "..\include\MultiBspline.h"
#include <iostream>
#include <vector>

class BsplineCurve : public MultiBspline
{
public:
    BsplineCurve(const std::vector<Bspline> &new_bspline, const std::vector<std::vector<double>> &new_controlPoints)
        : MultiBspline(new_bspline), controlPoints(new_controlPoints) {}

    ~BsplineCurve() {}

    Vertex<double> evaluateAtPoint(const double point);
    void plot2D(const int resolution, const std::string &filename) override;
    void plot3D(const int resolution, const std::vector<double> &zCoordinate, const std::string &filename);

    std::vector<std::vector<double>> &getControlPoints();

private:
    std::vector<std::vector<double>> controlPoints;
};

#include "..\src\BsplineCurve.cpp"

#endif