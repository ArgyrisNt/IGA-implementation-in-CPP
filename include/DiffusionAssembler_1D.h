#pragma once

#include <iostream>
#include "..\include\Matrix.h"
#include "..\include\Assembler_1D.h"

class DiffusionAssembler_1D : public Assembler_1D
{
public:
    // Constructor
    DiffusionAssembler_1D(double src, BoundCond &bcinfo, BsplineCurve &curve, double k, double delta)
        : Assembler_1D(src, bcinfo, curve), coef(k), Dt(delta) {}

    // Destructor
    virtual ~DiffusionAssembler_1D() {}

    // Member functions
    void assemble() override
    {
        calcStiff();
        calcMass();
        calcRhs();
        calcBound();
        sysMat = mass + stiff * (coef * Dt);
        for (int i = 0; i < rhs.size(); i++)
        {
            rhs[i] = rhs[i] * Dt;
        }
    }
    std::vector<double> nextStep(std::vector<double>);
    std::vector<double> applyInitCond(double (*func)(double));
    void applyBoundEllimination() override;
    void applyBoundMultipliers() override;
    void enforceBoundary(std::string&) override;

    // Member getter functions
    Matrix<double>& getMass() { return mass; }
    Matrix<double>& getSysMat() { return sysMat; }    

private:
    // Member local functions
    void calcMass();

    // Member variables
    Matrix<double> mass;
    Matrix<double> sysMat;
    double coef;
    double Dt; 
};