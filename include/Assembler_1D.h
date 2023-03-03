#pragma once

#include <iostream>
#include "..\include\Assembler.h"

class Assembler_1D : public Assembler
{
public:
    // Constructor
    Assembler_1D(double src, BoundCond &bcinfo, BsplineCurve &curve)
        : Assembler(src, bcinfo, curve.bspline_x), nOF(curve.bspline_x.getNOF()), ctrlPts(curve.ctrlPts) {}

    // Destructor
    virtual ~Assembler_1D() {}

    // Member functions
    void assemble() override
    {
        calcStiff();
        calcRhs();
        calcBound();
    }

    double calcJacobian(double, int, std::vector<double>&);
    
    // Member getter functions
    const int getNOF() const { return nOF; }
    std::vector<std::vector<double>> &getCtrlPts() { return ctrlPts; }

protected:
    // Member local functions 
    void calcStiff();
    void calcRhs();
    void calcBound();

    // Member variables
    const int nOF;
    std::vector<std::vector<double>> ctrlPts;
};