#pragma once

#include <iostream>
#include "..\include\Assembler.h"

class Assembler_1D : public Assembler
{
public:
    // Constructor
    Assembler_1D(double sourceFunction, BoundCond &boundaryConditions, BsplineCurve &curve)
        : Assembler(sourceFunction, boundaryConditions, curve.bspline_x), numberOfBasisFunctions(curve.bspline_x.getNumberOfBasisFunctions()), controlPoints(curve.controlPoints) {}

    // Destructor
    virtual ~Assembler_1D() {}

    // Member functions
    void assemble() override
    {
        computeStiffnessMatrix();
        computeRightHandSide();
        computeBoundary();
    }

    double Jacobian(double, int, std::vector<double>&);
    
    // Member getter functions
    const int getNumberOfBasisFunctions() const { return numberOfBasisFunctions; }
    std::vector<std::vector<double>> &getControlPoints() { return controlPoints; }

protected:
    // Member local functions 
    void computeStiffnessMatrix();
    void computeRightHandSide();
    void computeBoundary();

    // Member variables
    const int numberOfBasisFunctions;
    std::vector<std::vector<double>> controlPoints;
};