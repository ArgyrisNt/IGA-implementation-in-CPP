#ifndef H_ELEMENT
#define H_ELEMENT

#include <iostream>
#include <vector>
#include "..\include\TrimmingCurve.h"

class Element
{
public:
    Element(bool _isTrimmed, TrimmingCurve& _trimmingCurve) : isTrimmed(_isTrimmed), trimmingCurve(_trimmingCurve) {}

    void categorise();
    void computeTrimmedAndUntrimmedVertices();
    std::vector<std::vector<std::pair<double, double>>> divideInTriangles();
    std::vector<std::vector<std::pair<double, double>>> construct_3_triangles();
    std::vector<std::vector<std::pair<double, double>>> construct_2_triangles();
    std::vector<std::vector<std::pair<double, double>>> construct_1_triangle();

    bool isTrimmed;
    std::vector<std::pair<double, double>> untrimmedVertices;
    std::vector<std::pair<double, double>> trimmedVertices;
    std::vector<std::pair<double, double>> vertices;
    std::pair<double, double> centroid;
    TrimmingCurve trimmingCurve;
};

#include "..\src\Element.cpp"

#endif
