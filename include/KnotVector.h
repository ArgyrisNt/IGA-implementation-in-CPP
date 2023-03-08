#include <iostream>
#include <vector>

template<class T>
class KnotVector
{
public:
    KnotVector();
    KnotVector(int newDegree, std::vector<T>& newValues);
    KnotVector(const T start, const T end, int new_degree, int numberOfElements);

    KnotVector &operator=(const KnotVector &);
    T operator()(int position);

    int findSpanOfValue(const double value);
    void computeDistinctKnots(); // compute discrete knots
    void insert(int position, double value);

    int getSize() { return values.size(); }
    int getDegree() { return degree; }

    std::vector<double> distinctKnots;

private:
    std::vector<T> values;
    int degree;
};

#include "..\src\KnotVector.cpp"