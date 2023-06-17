#ifndef H_KNOTVECTOR
#define H_KNOTVECTOR

#include <iostream>
#include <vector>
#include "Utilities.h"

using namespace Utils;

template<class T>
class KnotVector
{
public:
    KnotVector() : degree(0) { computeDistinctKnots(); }
    KnotVector(const int newDegree, const std::vector<T> &newKnots, const std::vector<T> &newWeights);
    KnotVector(const T start, const T end, const int new_degree, const int numberOfElements, const std::vector<T> &newWeights);

    ~KnotVector() {}

    KnotVector &operator=(const KnotVector &);
    T operator()(int position);

    void computeDistinctKnots();
    int findSpanOfValue(const T value) const;
    void insert(const int position, const T value);
    std::vector<T> linspace(const int resolution) const;

    const int getSize() const 
    { return knots.size(); }

    const int getDegree() const 
    { return degree; }

    const std::vector<T> &getWeights() const 
    { return weights; }

    std::vector<T> &getDistinctKnots() 
    { return distinctKnots; }

    void setWeights(const std::vector<T>& new_weights);

private:
    std::vector<T> knots;
    std::vector<T> weights;
    std::vector<T> distinctKnots;
    int degree;
};

template <class T>
inline std::ostream &operator<<(std::ostream &os, KnotVector<T>& knotVector)
{
    for (int i = 0; i < knotVector.getSize(); ++i)
    {
        os << knotVector(i) << " ";
    }
    std::cout << std::endl;

    return os;
}

template <class T>
inline T KnotVector<T>::operator()(int position)
{
    assert(position >= 0);
    assert(position < knots.size());

    return knots[position];
}

#include "..\src\KnotVector.cpp"

#endif