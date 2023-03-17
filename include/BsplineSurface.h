#ifndef H_BSPLINESURFACE
#define H_BSPLINESURFACE

#include "..\include\MultiBspline.h"
#include "..\include\TrimmingCurve.h"
#include <iostream>
#include <vector>

class BsplineSurface : public MultiBspline
{
public:
    BsplineSurface(std::vector<Bspline> &&newBsplines, std::vector<std::vector<double>> &new_controlPoints, TrimmingCurve &_trimmingCurve)
        : MultiBspline(newBsplines), controlPoints(new_controlPoints), trimmingCurve(_trimmingCurve) {}

    ~BsplineSurface() {}

    BsplineSurface &operator=(const BsplineSurface &);

    Vertex<double> evaluateAtPoint(const Vertex<double> &&point);

    void plot2D(const int resolution, std::string filename) override;
    void plot3D(const int resolution, std::vector<double> &zCoordinate, std::string filename);

    void knotInsertion(KnotVector<double> &vector, std::vector<std::vector<double>> &points, const double newKnot);
    void uniformRefine_x();
    void uniformRefine_y();

    void setControlPoints(std::vector<std::vector<double>> &new_controlPoints);

    std::vector<std::vector<double>> &getControlPoints();

    TrimmingCurve trimmingCurve;

private:
    std::vector<std::vector<double>> XparametricCurvePoints(const int level);
    std::vector<std::vector<double>> YparametricCurvePoints(const int level);
    void refineParametricCurve(KnotVector<double> &vector, std::vector<std::vector<double>> &points);

    std::vector<std::vector<double>> controlPoints;
};

#include "..\src\BsplineSurface.cpp"

#endif