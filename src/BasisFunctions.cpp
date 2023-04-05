#include <iostream>
#include "..\include\BasisFunctions.h"


std::pair<std::vector<double>, std::vector<double>> BasisFunctions::evaluateAtPoint(const double value, bool allBasisFunctions)
{
    initializeBasisFunctions(value);

    for (int i = 1; i < knotVector.getDegree() + 1; ++i)
    {
        basisFunctionsOfDegree(i, value);
    }
    convertToNurbsFunctions();

    if (!allBasisFunctions)
    {
        computeActiveNurbsFunctions(value);
    }

    return std::make_pair(values, derivatives);
}

void BasisFunctions::basisFunctionsOfDegree(const int level, const double value)
{
    std::vector<double> oldValues = values;
    for (int j = 0; j < numberOfBasisFunctions; ++j)
    {
        double firstCoefficient = 0.0, secondCoefficient = 0.0;

        double denominator = knotVector(j + level) - knotVector(j);
        if (!almostEqual(denominator, 0.0))
            firstCoefficient = (value - knotVector(j)) / denominator;

        denominator = knotVector(j + level + 1) - knotVector(j + 1);
        if (!almostEqual(denominator, 0.0))
            secondCoefficient = (knotVector(j + level + 1) - value) / denominator;

        if (j == knotVector.getSize() - knotVector.getDegree() - 2)
            values[j] = firstCoefficient * oldValues[j];
        else
            values[j] = firstCoefficient * oldValues[j] + secondCoefficient * oldValues[j + 1];

        if (level == knotVector.getDegree())
        {
            denominator = knotVector(j + level) - knotVector(j);
            if (!almostEqual(denominator, 0.0))
                firstCoefficient = oldValues[j] / (knotVector(j + level) - knotVector(j));
            denominator = knotVector(j + level + 1) - knotVector(j + 1);
            if (!almostEqual(denominator, 0.0))
                secondCoefficient = oldValues[j + 1] / (knotVector(j + level + 1) - knotVector(j + 1));
            derivatives[j] = knotVector.getDegree() * (firstCoefficient - secondCoefficient);
        }
    }
}

void BasisFunctions::convertToNurbsFunctions()
{
    double sumValues = 0.0, sumDerivatives = 0.0;
    for (int i = 0; i < values.size(); ++i)
    {
        sumValues += values[i] * knotVector.getWeights()[i];
        sumDerivatives += derivatives[i] * knotVector.getWeights()[i];
    }

    std::vector<double> nurbs, nurbsDerivatives;
    for (int i = 0; i < values.size(); ++i)
    {
        double weightedBasis = values[i] * knotVector.getWeights()[i];
        double weightedDerivatives = derivatives[i] * knotVector.getWeights()[i];
        nurbs.push_back(weightedBasis / sumValues);
        nurbsDerivatives.push_back((weightedDerivatives * sumValues - weightedBasis * sumDerivatives) / (sumValues * sumValues));
    }

    values = nurbs;
    derivatives = nurbsDerivatives;
}

void BasisFunctions::computeActiveNurbsFunctions(const double value)
{
    // if value in [ui,ui+1), then active functions are Ni-p,...,Ni
    std::vector<double> activeValues;
    std::vector<double> activeDerivatives;
    int span = knotVector.findSpanOfValue(value);
    for (int i = span - knotVector.getDegree(); i <= span; ++i)
    {
        activeValues.push_back(values[i]);
        activeDerivatives.push_back(derivatives[i]);
    }

    values = activeValues;
    derivatives = activeDerivatives;
}

void BasisFunctions::initializeBasisFunctions(const double value)
{
    values = {};
    derivatives = {};
    numberOfBasisFunctions = knotVector.getSize() - knotVector.getDegree() - 1;
    std::vector<double> valuesOfBasisFunctions(numberOfBasisFunctions, 0.0);
    for (int j = 0; j < numberOfBasisFunctions; ++j)
    {
        bool isBetweenTheseTwoKnots = (value >= knotVector(j) && value < knotVector(j + 1));
        if (isBetweenTheseTwoKnots) valuesOfBasisFunctions[j] = 1.0;
    }

    bool isLastKnot = almostEqual(value, knotVector(knotVector.getSize() - 1));
    if (isLastKnot) valuesOfBasisFunctions[numberOfBasisFunctions - 1] = 1.0;

    values = valuesOfBasisFunctions;
    derivatives = valuesOfBasisFunctions;
}