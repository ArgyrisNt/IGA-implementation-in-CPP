#ifndef H_BSPLINECURVE
#define H_BSPLINECURVE

#include "BsplineEntity.h"
#include <iostream>
#include <vector>

class BsplineCurve : public BsplineEntity
{
public:
    BsplineCurve(const std::vector<Bspline> &new_bspline, const std::vector<Vertex<double>> &new_controlPoints)
        : BsplineEntity(MultiBspline(new_bspline), new_controlPoints) {}

    ~BsplineCurve() {}

    Vertex<double> evaluateAtPoint(const double point);
    void plot(const int resolution, const std::string &filename) override;
    void plotControlPoints(const std::string &filename);
    void plotVectorOnEntity(const int resolution, const std::vector<double> &zCoordinate, const std::string &filename) override;
    
    void uniformRefine_x() override;
};

#include "..\src\BsplineCurve.cpp"

#endif