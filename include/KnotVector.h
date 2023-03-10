#include <iostream>
#include <vector>

template<class T>
class KnotVector
{
public:
    KnotVector();
    KnotVector(int newDegree, std::vector<T>& newValues);
    KnotVector(const T start, const T end, int new_degree, int numberOfElements);

    ~KnotVector() {}

    KnotVector &operator=(const KnotVector &);
    T operator()(int position);

    int findSpanOfValue(const double value);
    void insert(int position, double value);

    int getSize() { return values.size(); }
    int getDegree() { return degree; }

private:
    std::vector<T> values;
    int degree;
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