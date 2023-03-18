#ifndef H_KNOTVECTOR
#define H_KNOTVECTOR

#include <iostream>
#include <vector>
#include "..\include\Utilities.h"

using namespace Utils;

template<class T>
class KnotVector
{
public:
    KnotVector();
    KnotVector(const int newDegree, const std::vector<T> &newKnots, const std::vector<T> &newWeights);
    KnotVector(const T start, const T end, const int new_degree, const int numberOfElements, const std::vector<T> &newWeights);

    ~KnotVector() {}

    KnotVector &operator=(const KnotVector &);
    T operator()(int position);

    void computeDistinctKnots();
    int findSpanOfValue(const T value) const;
    void insert(const int position, const T value);
    std::vector<T> linspace(const int resolution) const;

    int getSize();
    int getDegree();
    const std::vector<T> &getWeights() const;
    std::vector<T>& getDistinctKnots();

    void setWeights(const std::vector<T>& new_weights);

private:
    std::vector<T> knots;
    std::vector<T> weights;
    int degree;
    std::vector<T> distinctKnots;
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

#include "..\src\KnotVector.cpp"

#endif