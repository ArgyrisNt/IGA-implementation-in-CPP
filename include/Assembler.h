#pragma once

#include <iostream>
#include "..\include\Matrix.h"
#include "..\include\BoundCond.h"
#include "..\include\Bspline.h"

class Assembler
{
public:
    // Cunstructor
    Assembler(double newSourceFunction, BoundCond &newBoundaryConditions, Bspline &bspline)
        : sourceFunction(newSourceFunction), boundaryConditions(&newBoundaryConditions), bspline_x(&bspline) {}

    // Destructor
    virtual ~Assembler() {}

    // Member functions
    virtual void assemble() = 0;
    virtual void applyBoundaryEllimination();
    virtual void applyBoundaryMultipliers();
    virtual void enforceBoundaryConditions(std::string&);

    // Member getter functions
    Matrix<double> &getStiffnessMatrix() { return stiffnessMatrix; }
    std::vector<double> &getRightHandSide() { return rightHandSide; }
    Bspline& getBspline_x() { return *bspline_x; }
    BoundCond &getBoundaryConditions() { return *boundaryConditions; }
    std::string& getBoundaryMode() { return boundaryMode; }
    std::vector<std::pair<int, int>> &getBoundaryBasisFunctions() { return boundaryBasisFunctions; }

protected:
    // Member variables
    Matrix<double> stiffnessMatrix;
    std::vector<double> rightHandSide;
    Bspline* bspline_x;
    double sourceFunction;
    BoundCond* boundaryConditions;
    std::string boundaryMode;
    std::vector<std::pair<int, int>> boundaryBasisFunctions;
};