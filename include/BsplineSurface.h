#ifndef H_BSPLINESURFACE
#define H_BSPLINESURFACE

#include "..\include\Bspline.h"
#include <iostream>
#include <vector>

class BsplineSurface
{
public:
    BsplineSurface(const Bspline &new_bspline_x, const Bspline &new_bspline_y, std::vector<std::vector<double>> &new_controlPoints, 
                    TrimmingCurve& _trimmingCurve)
        : bspline_x(new_bspline_x), bspline_y(new_bspline_y), controlPoints(new_controlPoints), trimmingCurve(_trimmingCurve) {}

    ~BsplineSurface() {}

    BsplineSurface &operator=(const BsplineSurface &);

    Vertex<double> evaluateAtPoint(Vertex<double> &&point);

    void plot2D(int resolution, std::string filename);
    void plot3D(int resolution, std::vector<double>& zCoordinate, std::string filename);

    void knotInsertion(KnotVector<double> &vector, std::vector<std::vector<double>> &points, double newKnot);
    void uniformRefine_x();
    void uniformRefine_y();

    void setControlPoints(std::vector<std::vector<double>> &new_controlPoints);
    void setBspline_x(Bspline &new_bspline_x);
    void setBspline_y(Bspline &new_bspline_y);

    Bspline &getBspline_x();
    Bspline &getBspline_y();
    std::vector<std::vector<double>> &getControlPoints();

    TrimmingCurve trimmingCurve;

private:
    std::vector<std::vector<double>> XparametricCurvePoints(int level);
    std::vector<std::vector<double>> YparametricCurvePoints(int level);
    void refineParametricCurve(KnotVector<double> &vector, std::vector<std::vector<double>> &points);

    Bspline bspline_x;
    Bspline bspline_y;
    std::vector<std::vector<double>> controlPoints;
};

#include "..\src\BsplineSurface.cpp"

#endif