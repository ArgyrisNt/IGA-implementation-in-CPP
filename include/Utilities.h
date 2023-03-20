#ifndef H_UTILITIES
#define H_UTILITIES

#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include "math.h"
#include "Triangle.h"


namespace Utils
{

    struct CompareFirst
    {
        CompareFirst(int val) : val_(val) {}
        bool operator()(const std::pair<int, char> &elem) const
        {
            return val_ == elem.first;
        }

    private:
        int val_;
    };

    template <class T>
    std::vector<T> createTensorProduct(std::vector<T> &vec1, std::vector<T> &vec2);

    template <class T>
    void mapValuesToDomain(std::vector<T> &GaussPoints, const T left, const T right);

    template <class T>
    std::vector<std::pair<double, double>> GaussPointsAndWeightsQuad(int numberOfPoints, const double left, const double right);

    template <class T>
    std::vector<std::pair<double, double>> GaussPointsAndWeightsTriangle(double a, double b);

    bool almostEqual(const double a, const double b)
    { return fabs(a - b) < 1e-7; }

    template <class T>
    T norm(const std::vector<T> &v);

    template <class T>
    std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec);

} // end namespace Utilities

#include "..\src\Utilities.cpp"

#endif