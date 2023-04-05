#ifndef H_BSPLINEENTITY
#define H_BSPLINEENTITY

#include <iostream>
#include <vector>
#include "MultiBspline.h"

class BsplineEntity
{
public:
    BsplineEntity(const MultiBspline& _multiBspline, const std::vector<Vertex<double>>& _controlPoints)
        : multiBspline(_multiBspline), controlPoints(_controlPoints) {}

    virtual ~BsplineEntity() {}

    MultiBspline& getMultiBspline()
    { return multiBspline; }

    const std::vector<Vertex<double>> &getControlPoints() const 
    { return controlPoints; }

    void setControlPoints(const std::vector<Vertex<double>> &new_controlPoints) 
    { controlPoints = new_controlPoints; }

    void plotControlPoints(const std::string &filename);
    virtual void plot(const int, const std::string &) = 0;
    virtual void plotVectorOnEntity(const int, const std::vector<double> &, const std::string &) = 0;
    virtual void uniformRefine_x() = 0;

protected:
    MultiBspline multiBspline;
    std::vector<Vertex<double>> controlPoints;

    void knotInsertion(KnotVector<double> &vector, std::vector<Vertex<double>> &points, const double newKnot);
};

#include "..\src\BsplineEntity.cpp"

#endif