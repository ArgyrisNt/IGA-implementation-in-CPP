#ifndef H_BASISFUNCTIONS
#define H_BASISFUNCTIONS

#include <iostream>
#include <vector>
#include "KnotVector.h"

class BasisFunctions
{
public:
    BasisFunctions() {}
    BasisFunctions(KnotVector<double>& newKnotVector) : knotVector(newKnotVector) {}

    // all basis functions if all = true
    std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(const double value, bool all = false);
    void initializeBasisFunctions(const double value);
    void basisFunctionsOfDegree(const int level, const double value);
    void computeActiveBasisFunctions(const double value);

    int numberOfBasisFunctions;
    std::vector<double> values;
    std::vector<double> derivatives;
    KnotVector<double> knotVector;
};

#include "..\src\BasisFunctions.cpp"

#endif