#ifndef H_DIFFUSIONASSEMBLER
#define H_DIFFUSIONASSEMBLER

#include <iostream>
#include "..\include\Matrix.h"
#include "..\include\Assembler_1D.h"

class DiffusionAssembler_1D : public Assembler_1D
{
public:
    DiffusionAssembler_1D(const double newSourceFunction, const BoundCond &boundaryConditions, const BsplineCurve &curve, 
            const double newCoefficient, const double delta)
        : Assembler_1D(newSourceFunction, boundaryConditions, curve), coefficient(newCoefficient), Timestep(delta) {}

    virtual ~DiffusionAssembler_1D() {}

    void assemble() override;

    std::vector<double> nextStep(const std::vector<double> &);
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