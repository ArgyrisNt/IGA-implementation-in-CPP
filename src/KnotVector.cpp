//#include "..\include\KnotVector.h"
#include <iostream>
#include <cassert>
#include <algorithm>

template <class T>
KnotVector<T>::KnotVector()
{
    degree = 0;
    values = {};
}

template <class T>
KnotVector<T>::KnotVector(int newDegree, std::vector<T> &newValues)
{
    degree = newDegree;
    values = newValues;
}

template<class T>
KnotVector<T>::KnotVector(const T start, const T end, int newDegree, int numberOfElements)
{
    degree = newDegree;
    for (int i = 0; i < degree; i++)
    {
        values.push_back(start);
    }
    double delta = end / numberOfElements;
    for (int i = 0; i <= numberOfElements; i++)
    {
        values.push_back(i * delta);
    }
    for (int i = 0; i < degree; i++)
    {
        values.push_back(end);
    }
}

template <class T>
KnotVector<T> &KnotVector<T>::operator=(const KnotVector &oldKnotVector)
{
    values = oldKnotVector.values;
    degree = oldKnotVector.degree;
}

template <class T>
T KnotVector<T>::operator()(int position)
{
    assert(position >= 0);
    assert(position < values.size());
    return values[position];
}

template <class T>
int KnotVector<T>::findSpanOfValue(const double value)
{
    bool isLastKnot = (value == values[values.size() - 1]);
    if (isLastKnot) return (values.size() - degree - 2);
    for (int index = 0; index < values.size() - 1; index++)
    {
        if (value >= values[index] && value < values[index + 1]) return index;
    }
    std::cout << "Error: Value " << value << " is not in vector." << std::endl;
    throw std::invalid_argument("Value does not appear in vector");
}

template <class T>
void KnotVector<T>::insert(int position, double value)
{
    values.insert(values.begin() + position, value);
}
