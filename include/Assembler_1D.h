#ifndef H_ASSEMBLER_1D
#define H_ASSEMBLER_1D

#include <iostream>
#include "Assembler.h"
#include "BsplineCurve.h"

class Assembler_1D : public Assembler<BsplineCurve>
{
public:
    Assembler_1D(const double sourceFunction, const BoundCond &boundaryConditions, const BsplineCurve &_curve)
        : Assembler<BsplineCurve>(sourceFunction, boundaryConditions, _curve) {}

    ~Assembler_1D() {}

    const int getNumberOfBasisFunctions()
    { return getBspline_x().getNumberOfBasisFunctions(); }

    const std::vector<Vertex<double>> &getControlPoints() const
    { return bsplineEntity->getControlPoints();  }

    Bspline &getBspline_x()
    { return bsplineEntity->getMultiBspline().getBspline(0); }

protected:
    Matrix<double> Jacobian(double, const std::vector<double> &);

    void computeBoundary() override;
    void computeStiffnessMatrixAndRightHandSide() override;
    double computeStiffnessIntegral(int basisFunction, int trialFunction);
    double computeRightHandSideIntegral(int basisFunction);
};

#include "..\src\Assembler_1D.cpp"

#endif