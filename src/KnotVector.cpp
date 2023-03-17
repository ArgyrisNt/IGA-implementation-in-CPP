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
std::vector<T>& KnotVector<T>::getWeights()
{
    return weights;
}

template <class T>
std::vector<T>& KnotVector<T>::getDistinctKnots()
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
KnotVector<T>::KnotVector(const int newDegree, std::vector<T> &newKnots, std::vector<T> &newWeights)
{
    degree = newDegree;
    knots = newKnots;
    weights = newWeights;
    computeDistinctKnots();
}

template <class T>
KnotVector<T>::KnotVector(const T start, const T end, const int newDegree, const int numberOfElements, std::vector<T> &newWeights)
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
        if (!almostEqual(currentValue, previousValue)) distinctKnots.push_back(knots[i]);
        previousValue = knots[i];
    }
}

template <class T>
int KnotVector<T>::findSpanOfValue(const T value)
{
    bool isLastKnot = almostEqual(value, knots[knots.size() - 1]);
    if (isLastKnot) return (knots.size() - degree - 2);
    for (int index = 0; index < knots.size() - 1; index++)
    {
        if (value >= knots[index] && value < knots[index + 1]) return index;
    }
    std::cout << "Error: Value " << value << " is not in vector." << std::endl;
    throw std::invalid_argument("Value does not appear in vector");
}

template <class T>
void KnotVector<T>::insert(const int position, const T value)
{
    assert(position >= 0 && position <= knots.size());
    knots.insert(knots.begin() + position, value);
}

template <class T>
std::vector<T> KnotVector<T>::linspace(const int resolution)
{
    assert(resolution != 0);
    std::vector<T> steps;
    T start = knots[0];
    T end = knots[knots.size() - 1];
    for (int i = (int)(start); i <= resolution; i++)
    {
        T i_step = start + (T)(i) * ((end - start) / ((T)(resolution)));
        steps.push_back(i_step);
    }
    return steps;
}

template <class T>
void KnotVector<T>::setWeights(std::vector<T>& new_weights)
{
    for (auto weight : new_weights)
    {
        assert(weight >= 0.0 && weight <= 1.0);
    }

    weights = new_weights;
}