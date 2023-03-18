#ifndef H_BSPLINESURFACE
#define H_BSPLINESURFACE

#include "..\include\MultiBspline.h"
#include "..\include\TrimmingCurve.h"
#include <iostream>
#include <vector>

class BsplineSurface : public MultiBspline
{
public:
    BsplineSurface(const std::vector<Bspline> &newBsplines, const std::vector<std::vector<double>> &new_controlPoints, const TrimmingCurve &_trimmingCurve)
        : MultiBspline(newBsplines), controlPoints(new_controlPoints), trimmingCurve(_trimmingCurve) {}

    ~BsplineSurface() {}

    BsplineSurface &operator=(const BsplineSurface &);

    Vertex<double> evaluateAtPoint(const Vertex<double> &&point);

    void plot2D(const int resolution, const std::string &filename) override;
    void plot3D(const int resolution, const std::vector<double> &zCoordinate, const std::string &filename);

    void knotInsertion(KnotVector<double> &vector, std::vector<std::vector<double>> &points, const double newKnot);
    void uniformRefine_x();
    void uniformRefine_y();

    void setControlPoints(const std::vector<std::vector<double>> &new_controlPoints);

    std::vector<std::vector<double>> &getControlPoints();
    TrimmingCurve &getTrimmingCurve();

private:
    std::vector<std::vector<double>> XparametricCurvePoints(const int level);
    std::vector<std::vector<double>> YparametricCurvePoints(const int level);
    void refineParametricCurve(KnotVector<double> &vector, std::vector<std::vector<double>> &points);

    std::vector<std::vector<double>> controlPoints;
    TrimmingCurve trimmingCurve;
};

#include "..\src\BsplineSurface.cpp"

#endif