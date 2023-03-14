#include "..\include\KnotVector.h"
#include <iostream>
#include <cassert>
#include <algorithm>

template <class T>
int KnotVector<T>::getSize()
{
    return knots.size();
}

template <class T>
int KnotVector<T>::getDegree()
{
    return degree;
}

template <class T>
std::vector<T> KnotVector<T>::getWeights()
{
    return weights;
}

template <class T>
std::vector<T> KnotVector<T>::getDistinctKnots()
{
    return distinctKnots;
}

template <class T>
KnotVector<T>::KnotVector()
{
    degree = 0;
    knots = {};
    weights = {};
    computeDistinctKnots();
}

template <class T>
KnotVector<T>::KnotVector(int newDegree, std::vector<T> &newKnots, std::vector<T> &newWeights)
{
    degree = newDegree;
    knots = newKnots;
    weights = newWeights;
    computeDistinctKnots();
}

template <class T>
KnotVector<T>::KnotVector(const T start, const T end, int newDegree, int numberOfElements, std::vector<T> &newWeights)
{
    degree = newDegree;
    for (int i = 0; i < degree; i++)
    {
        knots.push_back(start);
    }
    double delta = end / numberOfElements;
    for (int i = 0; i <= numberOfElements; i++)
    {
        knots.push_back(i * delta);
    }
    for (int i = 0; i < degree; i++)
    {
        knots.push_back(end);
    }
    weights = newWeights;
    computeDistinctKnots();
}

template <class T>
KnotVector<T> &KnotVector<T>::operator=(const KnotVector &oldKnotVector)
{
    knots = oldKnotVector.knots;
    degree = oldKnotVector.degree;
    weights = oldKnotVector.weights;
    distinctKnots = oldKnotVector.distinctKnots;
}

template <class T>
T KnotVector<T>::operator()(int position)
{
    assert(position >= 0);
    assert(position < knots.size());

    return knots[position];
}

template <class T>
void KnotVector<T>::computeDistinctKnots()
{
    distinctKnots = {};
    double currentValue, previousValue = -100.0;
    for (int i = degree; i < getSize() - degree; i++)
    {
        currentValue = knots[i];
        if (currentValue != previousValue)
            distinctKnots.push_back(knots[i]);
        previousValue = knots[i];
    }
}

template <class T>
int KnotVector<T>::findSpanOfValue(const double value)
{
    bool isLastKnot = (value == knots[knots.size() - 1]);
    if (isLastKnot) return (knots.size() - degree - 2);
    for (int index = 0; index < knots.size() - 1; index++)
    {
        if (value >= knots[index] && value < knots[index + 1]) return index;
    }
    std::cout << "Error: Value " << value << " is not in vector." << std::endl;
    throw std::invalid_argument("Value does not appear in vector");
}

template <class T>
void KnotVector<T>::insert(int position, double value)
{
    knots.insert(knots.begin() + position, value);
}

template <class T>
void KnotVector<T>::setWeights(std::vector<double>& new_weights)
{
    for (auto weight : new_weights)
    {
        if (weight < 0 || weight > 1)
        {
            std::cout << "Invalid weights. Valid weights are within [0,1]." << std::endl;
            throw std::invalid_argument("Invalid weights");
            break;
        }
    }

    weights = new_weights;
}