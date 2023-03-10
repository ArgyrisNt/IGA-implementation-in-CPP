#pragma once

#include <iostream>
#include "..\include\Assembler.h"

class Assembler_1D : public Assembler
{
public:
    Assembler_1D(double sourceFunction, BoundCond &boundaryConditions, BsplineCurve &curve)
        : Assembler(sourceFunction, boundaryConditions, curve.getBspline_x()), numberOfBasisFunctions(curve.getBspline_x().getNumberOfBasisFunctions()), controlPoints(curve.getControlPoints()) {}

    virtual ~Assembler_1D() {}

    void assemble() override;

    Matrix<double> Jacobian(double, std::vector<double> &);

    const int getNumberOfBasisFunctions() const { return numberOfBasisFunctions; }
    std::vector<std::vector<double>> &getControlPoints() { return controlPoints; }

protected:
    void computeBoundary();
    void computeStiffnessMatrixAndRightHandSide();
    double computeStiffnessIntegral(int element, int basisFunction, int trialFunction);
    double computeRightHandSideIntegral(int element, int basisFunction);

    const int numberOfBasisFunctions;
    std::vector<std::vector<double>> controlPoints;
};