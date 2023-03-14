#ifndef H_KNOTVECTOR
#define H_KNOTVECTOR

#include <iostream>
#include <vector>

template<class T>
class KnotVector
{
public:
    KnotVector();
    KnotVector(int newDegree, std::vector<T> &newKnots, std::vector<T> &newWeights);
    KnotVector(const T start, const T end, int new_degree, int numberOfElements, std::vector<T> &newWeights);

    ~KnotVector() {}

    KnotVector &operator=(const KnotVector &);
    T operator()(int position);

    void computeDistinctKnots();
    int findSpanOfValue(const double value);
    void insert(int position, double value);

    int getSize();
    int getDegree();
    std::vector<T> getWeights();
    std::vector<T> getDistinctKnots();

    void setWeights(std::vector<double>& new_weights);

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