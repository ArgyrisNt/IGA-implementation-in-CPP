#ifndef H_ELEMENT
#define H_ELEMENT

#include <iostream>
#include <vector>
#include "..\include\TrimmingCurve.h"

class Element
{
public:
    Element(bool _isTrimmed, TrimmingCurve& _trimmingCurve) : isTrimmed(_isTrimmed), trimmingCurve(_trimmingCurve) {}

    ~Element() {}

    void categorise();
    void computeTrimmedAndUntrimmedVertices();
    std::vector<Triangle<double>> divideInTriangles();
    std::vector<Triangle<double>> construct_3_triangles();
    std::vector<Triangle<double>> construct_2_triangles();
    std::vector<Triangle<double>> construct_1_triangle();

    void setVertices(std::vector<Vertex<double>>& _vertices) { vertices =_vertices; }

    bool isTrimmed;

private:
    Vertex<double> centroid;
    TrimmingCurve trimmingCurve;
    std::vector<Vertex<double>> vertices;
    std::vector<Vertex<double>> trimmedVertices;
    std::vector<Vertex<double>> untrimmedVertices;
};

#include "..\src\Element.cpp"

#endif
