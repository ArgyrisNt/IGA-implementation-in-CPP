#ifndef H_ASSEMBLER
#define H_ASSEMBLER

#include <iostream>
#include "Matrix.h"
#include "BoundCond.h"
#include "Bspline.h"

template <class T>
class Assembler
{
public:
    Assembler() {}
    Assembler(const double newSourceFunction, const BoundCond &newBoundaryConditions, const T &bspline)
        : sourceFunction(newSourceFunction), boundaryConditions(newBoundaryConditions), bsplineEntity(std::make_shared<T>(bspline)) {}

    virtual ~Assembler() {}

    virtual void assemble();
    void enforceBoundaryConditions(const std::string &);

    int spanOfValueInKnotVector(const double value, const int dim)
    { return bsplineEntity->getMultiBspline().findSpanOfValue(value, dim); }

    Matrix<double> &getSystemMatrix() 
    { return systemMatrix; }

    std::vector<double> &getRightHandSide() 
    { return rightHandSide; }

    std::vector<double> &getDistinctKnots(const int dim);
    const std::string &getBoundarymode() const 
    { return boundaryMode; }

    const BoundCond &getBoundaryConditions() const 
    { return boundaryConditions; }

    const std::vector<std::pair<int, int>> &getBoundaryBasisFunctions() const 
    { return boundaryBasisFunctions; }

protected:
    virtual void computeBoundary() = 0;
    virtual void computeStiffnessMatrixAndRightHandSide() = 0;

    void applyBoundaryEllimination();
    void applyBoundaryMultipliers();
    double addBoundaryValueToRhs(const int position);

    Matrix<double> stiffnessMatrix;
    Matrix<double> systemMatrix;
    std::vector<double> rightHandSide;

    const double sourceFunction;
    std::vector<std::pair<double, double>> XGaussPointsAndWeights;

    std::string boundaryMode;
    BoundCond boundaryConditions;
    std::vector<std::pair<int, int>> boundaryBasisFunctions;

    std::shared_ptr<T> bsplineEntity;
};

template <class T>
inline std::vector<double> &Assembler<T>::getDistinctKnots(const int dim)
{
    assert(dim >= 0 && dim < bsplineEntity->getMultiBspline().getDimension());
    return bsplineEntity->getMultiBspline().getKnotvector(dim).getDistinctKnots();
}

#include "..\src\Assembler.cpp"

#endif