//#include "..\include\KnotVector.h"
#include <iostream>
#include <cassert>
#include <algorithm>

template <class T>
KnotVector<T>::KnotVector()
{
    degree = 0;
    values = {};
    computeDistinctKnots();
}

template <class T>
KnotVector<T>::KnotVector(int newDegree, std::vector<T> &newValues)
{
    degree = newDegree;
    values = newValues;
    computeDistinctKnots();
}

template<class T>
KnotVector<T>::KnotVector(const T start, const T end, int new_degree, int numberOfElements)
{
    degree = new_degree;
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
    computeDistinctKnots();
}

template <class T>
KnotVector<T> &KnotVector<T>::operator=(const KnotVector &oldKnotVector)
{
    values = oldKnotVector.values;
    degree = oldKnotVector.degree;
    distinctKnots = oldKnotVector.distinctKnots;
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
    int index = -1;
    if (value == values[values.size() - 1])
    {
        return (values.size() - degree - 2);
    }
    for (int i = 0; i < values.size() - 1; i++)
    {
        double start = values[i];
        double end = values[i + 1];
        if (value >= start && value < end)
        {
            index = i;
            break;
        }
    }
    if (index == -1)
    {
        std::cout << "Error: Value " << value << " is not in vector." << std::endl;
        throw std::invalid_argument("Value does not appear in vector");
    }

    return index;
}

template <class T>
void KnotVector<T>::computeDistinctKnots()
{
    distinctKnots = {};
    for (int i = degree; i < getSize() - degree; i++)
    {
        if (!std::count(distinctKnots.begin(), distinctKnots.end(), values[i]))
        {
            distinctKnots.push_back(values[i]);
        }
    }
}

template <class T>
void KnotVector<T>::insert(int position, double value)
{
    values.insert(values.begin() + position, value);
}
