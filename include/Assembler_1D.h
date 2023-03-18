#ifndef H_ASSEMBLER_1D
#define H_ASSEMBLER_1D

#include <iostream>
#include "..\include\Assembler.h"
#include "..\include\BsplineCurve.h"

class Assembler_1D : public Assembler<BsplineCurve>
{
public:
    Assembler_1D(const double sourceFunction, const BoundCond &boundaryConditions, const BsplineCurve &_curve)
        : Assembler<BsplineCurve>(sourceFunction, boundaryConditions, _curve) {}

    ~Assembler_1D() {}

    void assemble() override;

    Matrix<double> Jacobian(double, const std::vector<double> &);

    const int getNumberOfBasisFunctions();
    std::vector<std::vector<double>> &getControlPoints();
    Bspline &getBspline_x();

protected:
    void computeBoundary();
    void computeStiffnessMatrixAndRightHandSide();
    double computeStiffnessIntegral(int element, int basisFunction, int trialFunction);
    double computeRightHandSideIntegral(int element, int basisFunction);
};

#include "..\src\Assembler_1D.cpp"

#endif