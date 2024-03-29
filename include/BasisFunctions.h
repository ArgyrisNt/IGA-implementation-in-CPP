#ifndef H_BASISFUNCTIONS
#define H_BASISFUNCTIONS

#include <iostream>
#include <vector>
#include <memory>
#include "KnotVector.h"

class BasisFunctions
{
public:
    BasisFunctions() {}
    BasisFunctions(const KnotVector<double>& newKnotVector) : knotVector(newKnotVector) {}

    // all basis functions if all = true
    std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(const double value, bool all = false);

private:
    void initializeBasisFunctions(const double value);
    void basisFunctionsOfDegree(const int level, const double value);
    void computeActiveNurbsFunctions(const double value);
    void convertToNurbsFunctions();

    int numberOfBasisFunctions;
    std::vector<double> values;
    std::vector<double> derivatives;
    KnotVector<double> knotVector;
};

#include "..\src\BasisFunctions.cpp"

#endif