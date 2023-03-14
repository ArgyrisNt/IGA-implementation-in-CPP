#ifndef H_BASISFUNCTIONS
#define H_BASISFUNCTIONS

#include <iostream>
#include <vector>
#include "KnotVector.h"

class BasisFunctions
{
public:
    // all basis functions if all = true
    std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(const double value, bool all = false); 
    void initializeBasisFunctions(double value);
    void basisFunctionsOfDegree(int level, double value);
    void computeActiveBasisFunctions(double value);

    void setKnotVector(KnotVector<double>& newKnotVector);

    int numberOfBasisFunctions;
    std::vector<double> values;
    std::vector<double> derivatives;
    KnotVector<double> knotVector;
};

#include "..\src\BasisFunctions.cpp"

#endif