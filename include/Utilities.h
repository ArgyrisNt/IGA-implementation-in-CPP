#ifndef H_UTILITIES
#define H_UTILITIES

#include <iostream>
#include "math.h"

template <class T>
class Vertex
{
public:
    Vertex() {}
    Vertex(T _x, T _y) : x(_x), y(_y) {}

    void set(T _x, T _y) { x = _x; y = _y; }

    T x;
    T y;
};

template <class T>
bool almostEqual(const Vertex<T> &vertex1, const Vertex<T> &vertex2)
{
    return (fabs(vertex1.x - vertex2.x) < 1e-7 && fabs(vertex1.y - vertex2.y) < 1e-7);
}

template <class T>
class Triangle
{
public:
    Triangle(Vertex<T> _vertex1, Vertex<T> _vertex2, Vertex<T> _vertex3) 
        : vertex1(_vertex1), vertex2(_vertex2), vertex3(_vertex3) {}
    
    ~Triangle() {}

    Vertex<T> vertex1;
    Vertex<T> vertex2;
    Vertex<T> vertex3;
};

#endif