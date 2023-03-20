#ifndef H_ELEMENT
#define H_ELEMENT

#include <iostream>
#include <vector>
#include "TrimmingCurve.h"

class Element
{
public:
    Element(const TrimmingCurve& _trimmingCurve) : trimmingCurve(_trimmingCurve) {}

    ~Element() {}

    void categorise();
    std::vector<Triangle<double>> divideInTriangles();
    std::vector<Triangle<double>> construct_3_triangles();
    std::vector<Triangle<double>> construct_2_triangles();
    std::vector<Triangle<double>> construct_1_triangle();

    void setVertices(std::vector<Vertex<double>>& _vertices) 
    { vertices =_vertices; }

    bool isTrimmed() 
    { return is_Trimmed; }

private:
    void computeTrimmedAndUntrimmedVertices();

    Vertex<double> centroid;
    TrimmingCurve trimmingCurve;
    std::vector<Vertex<double>> vertices;
    std::vector<Vertex<double>> trimmedVertices;
    std::vector<Vertex<double>> untrimmedVertices;
    bool is_Trimmed = false;
};

#include "..\src\Element.cpp"

#endif
