#ifndef H_BSPLINESURFACE
#define H_BSPLINESURFACE

#include "..\include\Bspline_2D.h"
#include "..\include\TrimmingCurve.h"
#include <iostream>
#include <vector>

class BsplineSurface : public Bspline_2D
{
public:
    BsplineSurface(const Bspline &new_bspline_x, const Bspline &new_bspline_y, std::vector<std::vector<double>> &new_controlPoints, 
                    TrimmingCurve& _trimmingCurve)
        : Bspline_2D(new_bspline_x, new_bspline_y), controlPoints(new_controlPoints), trimmingCurve(_trimmingCurve) {}

    ~BsplineSurface() {}

    BsplineSurface &operator=(const BsplineSurface &);

    Vertex<double> evaluateAtPoint(Vertex<double> &&point);

    void plot2D(int resolution, std::string filename);
    void plot3D(int resolution, std::vector<double>& zCoordinate, std::string filename);

    void knotInsertion(KnotVector<double> &vector, std::vector<std::vector<double>> &points, double newKnot);
    void uniformRefine_x();
    void uniformRefine_y();

    void setControlPoints(std::vector<std::vector<double>> &new_controlPoints);

    std::vector<std::vector<double>> &getControlPoints();

    TrimmingCurve trimmingCurve;

private:
    std::vector<std::vector<double>> XparametricCurvePoints(int level);
    std::vector<std::vector<double>> YparametricCurvePoints(int level);
    void refineParametricCurve(KnotVector<double> &vector, std::vector<std::vector<double>> &points);

    std::vector<std::vector<double>> controlPoints;
};

#include "..\src\BsplineSurface.cpp"

#endif