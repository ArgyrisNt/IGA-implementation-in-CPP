#pragma once

#include <iostream>
#include "..\include\Matrix.h"
#include "..\include\Assembler_1D.h"

class DiffusionAssembler_1D : public Assembler_1D
{
public:
    // Constructor
    DiffusionAssembler_1D(double newSourceFunction, BoundCond &boundaryConditions, BsplineCurve &curve, double newCoefficient, double delta)
        : Assembler_1D(newSourceFunction, boundaryConditions, curve), coefficient(newCoefficient), Timestep(delta) {}

    // Destructor
    virtual ~DiffusionAssembler_1D() {}

    // Member functions
    void assemble() override
    {
        XcomputeDistinctKnots();
        computeStiffnessMatrixAndRightHandSide();
        computeMassMatrix();
        computeBoundary();
        systemMatrix = massMatrix + stiffnessMatrix * (coefficient * Timestep);
        for (int i = 0; i < rightHandSide.size(); i++)
        {
            rightHandSide[i] = rightHandSide[i] * Timestep;
        }
    }
    std::vector<double> nextStep(std::vector<double>);
    std::vector<double> applyInitialCondition(double (*func)(double));
    void applyBoundaryEllimination() override;
    void applyBoundaryMultipliers() override;
    void enforceBoundaryConditions(std::string&) override;

    // Member getter functions
    Matrix<double> &getMassMatrix() { return massMatrix; }
    Matrix<double> &getSystemMatrix() { return systemMatrix; }

private:
    // Member local functions
    void computeMassMatrix();
    double computeMassIntegral(int element, int basisFunction, int trialFunction);


    // Member variables
    Matrix<double> massMatrix;
    Matrix<double> systemMatrix;
    double coefficient;
    double Timestep; 
};