#ifndef H_BSPLINESURFACE
#define H_BSPLINESURFACE

#include "BsplineEntity.h"
#include "TrimmingCurve.h"
#include <iostream>
#include <vector>

class BsplineSurface : public BsplineEntity
{
public:
    BsplineSurface(const std::vector<Bspline> &newBsplines, const std::vector<Vertex<double>> &new_controlPoints, const TrimmingCurve &_trimmingCurve)
        : BsplineEntity(MultiBspline(newBsplines), new_controlPoints), trimmingCurve(_trimmingCurve) {}

    ~BsplineSurface() {}

    BsplineSurface &operator=(const BsplineSurface &);

    Vertex<double> evaluateAtPoint(const Vertex<double> &point);
    void plot(const int resolution, const std::string &filename) override;
    void plotVectorOnEntity(const int resolution, const std::vector<double> &zCoordinate, const std::string &filename) override;

    void uniformRefine_x() override;
    void uniformRefine_y();

    const TrimmingCurve &getTrimmingCurve() const 
    { return trimmingCurve; }

private:
    std::vector<Vertex<double>> XparametricCurvePoints(const int level);
    std::vector<Vertex<double>> YparametricCurvePoints(const int level);
    void refineParametricCurve(KnotVector<double> &vector, std::vector<Vertex<double>> &points);

    TrimmingCurve trimmingCurve;
};

#include "..\src\BsplineSurface.cpp"

#endif