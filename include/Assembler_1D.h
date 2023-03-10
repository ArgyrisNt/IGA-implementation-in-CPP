#pragma once

#include <iostream>
#include "..\include\Assembler.h"

class Assembler_1D : public Assembler
{
public:
    // Constructor
    Assembler_1D(double sourceFunction, BoundCond &boundaryConditions, BsplineCurve &curve)
        : Assembler(sourceFunction, boundaryConditions, curve.getBspline_x()), numberOfBasisFunctions(curve.getBspline_x().getNumberOfBasisFunctions()), controlPoints(curve.getControlPoints()) {}

    // Destructor
    virtual ~Assembler_1D() {}

    // Member functions
    void assemble() override
    {
        computeStiffnessMatrixAndRightHandSide();
        computeBoundary();
    }

    Matrix<double> Jacobian(double, std::vector<double> &);

    // Member getter functions
    const int getNumberOfBasisFunctions() const { return numberOfBasisFunctions; }
    std::vector<std::vector<double>> &getControlPoints() { return controlPoints; }

protected:
    // Member local functions
    void computeStiffnessMatrixAndRightHandSide();
    void computeBoundary();
    double computeStiffnessIntegral(int element, int basisFunction, int trialFunction);
    double computeRightHandSideIntegral(int element, int basisFunction);

    // Member variables
    const int numberOfBasisFunctions;
    std::vector<std::vector<double>> controlPoints;
};