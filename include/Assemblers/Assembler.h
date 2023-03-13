#ifndef H_ASSEMBLER
#define H_ASSEMBLER

#include <iostream>
#include "..\include\Matrix.h"
#include "..\include\BoundCond.h"
#include "..\include\Bspline.h"

class Assembler
{
public:
    Assembler() {}
    Assembler(double newSourceFunction, Bspline &bspline) : sourceFunction(newSourceFunction), bspline_x(&bspline) {}

    virtual ~Assembler() {}

    virtual void assemble() = 0;

    Matrix<double> &getStiffnessMatrix();
    Matrix<double> &getSystemMatrix();
    std::vector<double> &getRightHandSide();
    Bspline &getBspline_x();
    double getDistinctKnotX(int position);
    std::vector<double> getDistinctKnotsX();

    int XspanOfValueInKnotVector(double value);

protected:
    std::vector<std::pair<double, double>> XGaussPointsAndWeights;

    Matrix<double> stiffnessMatrix;
    Matrix<double> systemMatrix;
    std::vector<double> rightHandSide;
    Bspline *bspline_x;
    double sourceFunction;
};

class AssemblerBoundary : public Assembler
{
public:
    AssemblerBoundary(double newSourceFunction, BoundCond &newBoundaryConditions, Bspline &bspline) 
        : boundaryConditions(&newBoundaryConditions), Assembler(newSourceFunction, bspline) {}

    void applyBoundaryEllimination();
    void applyBoundaryMultipliers();
    void enforceBoundaryConditions(std::string &);
    double addBoundaryValueToRhs(int position);

    BoundCond *boundaryConditions;
    std::string boundaryMode;
    std::vector<std::pair<int, int>> boundaryBasisFunctions;
};

#include "..\src\Assembler.cpp"

#endif