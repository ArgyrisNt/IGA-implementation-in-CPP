#ifndef H_DIFFUSIONASSEMBLER
#define H_DIFFUSIONASSEMBLER

#include <iostream>
#include "..\include\Matrix.h"
#include "..\include\Assembler_1D.h"

class DiffusionAssembler_1D : public Assembler_1D
{
public:
    DiffusionAssembler_1D(double newSourceFunction, BoundCond &boundaryConditions, BsplineCurve &curve, 
            double newCoefficient, double delta)
        : Assembler_1D(newSourceFunction, boundaryConditions, curve), coefficient(newCoefficient), Timestep(delta) {}

    virtual ~DiffusionAssembler_1D() {}

    void assemble() override;

    std::vector<double> nextStep(std::vector<double>&);
    std::vector<double> applyInitialCondition(double (*func)(double));

    Matrix<double> &getMassMatrix();

private:
    void computeMassMatrix();
    double computeMassIntegral(int element, int basisFunction, int trialFunction);

    Matrix<double> massMatrix;
    double coefficient;
    double Timestep; 
};

#include "..\src\DiffusionAssembler_1D.cpp"

#endif