#ifndef H_ASSEMBLER
#define H_ASSEMBLER

#include <iostream>
#include "..\include\Matrix.h"
#include "..\include\BoundCond.h"
#include "..\include\Bspline.h"

template <class T>
class Assembler
{
public:
    Assembler() {}
    Assembler(const double newSourceFunction, const BoundCond &newBoundaryConditions, const T &bspline)
        : sourceFunction(newSourceFunction), bsplineEntity(bspline), boundaryConditions(newBoundaryConditions) {}

    virtual ~Assembler() {}

    virtual void assemble() = 0;

    Matrix<double> &getStiffnessMatrix();
    Matrix<double> &getSystemMatrix();
    std::vector<double> &getRightHandSide();
    std::vector<double> &getDistinctKnots(const int dim);

    int spanOfValueInKnotVector(const double value, const int dim);

    void applyBoundaryEllimination();
    void applyBoundaryMultipliers();
    void enforceBoundaryConditions(const std::string &);
    double addBoundaryValueToRhs(const int position);

    std::string &getBoundarymode() { return boundaryMode; }
    BoundCond &getBoundaryConditions() { return boundaryConditions; }
    std::vector<std::pair<int, int>> &getBoundaryBasisFunctions() { return boundaryBasisFunctions; }

protected:
    std::vector<std::pair<double, double>> XGaussPointsAndWeights;

    Matrix<double> stiffnessMatrix;
    Matrix<double> systemMatrix;
    std::vector<double> rightHandSide;
    T bsplineEntity;
    const double sourceFunction;
    std::string boundaryMode;
    BoundCond boundaryConditions;
    std::vector<std::pair<int, int>> boundaryBasisFunctions;
};

#include "..\src\Assembler.cpp"

#endif