#ifndef H_KNOTVECTOR
#define H_KNOTVECTOR

#include <iostream>
#include <vector>

template<class T>
class KnotVector
{
public:
    KnotVector();
    KnotVector(const int newDegree, std::vector<T> &newKnots, std::vector<T> &newWeights);
    KnotVector(const T start, const T end, const int new_degree, const int numberOfElements, std::vector<T> &newWeights);

    ~KnotVector() {}

    KnotVector &operator=(const KnotVector &);
    T operator()(int position);

    void computeDistinctKnots();
    int findSpanOfValue(const T value);
    void insert(const int position, const T value);
    std::vector<T> linspace(const int resolution);

    int getSize();
    int getDegree();
    std::vector<T>& getWeights();
    std::vector<T>& getDistinctKnots();

    void setWeights(std::vector<T>& new_weights);

private:
    std::vector<T> knots;
    std::vector<T> weights;
    int degree;
    std::vector<T> distinctKnots;
};

template <class T>
inline std::ostream &operator<<(std::ostream &os, KnotVector<T>& knotVector)
{
    for (int i = 0; i < knotVector.getSize(); i++)
    {
        os << knotVector(i) << " ";
    }
    std::cout << std::endl;

    return os;
}

#include "..\src\KnotVector.cpp"

#endif