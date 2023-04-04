#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>

template <class T>
class Vertex
{
public:
    Vertex() {}
    Vertex(T _x, T _y) : x(_x), y(_y) {}

    void set(T _x, T _y) { x = _x; y = _y; }

    bool operator==(Vertex& vertex2)
    { return (fabs(x - vertex2.x) < 1e-7 && fabs(y - vertex2.y) < 1e-7); }

    bool operator!=(Vertex& vertex2)
    { return !((*this) == vertex2); }


    T x;
    T y;
};

template <class T>
Vertex<T> operator+(const Vertex<T>& v1, const Vertex<T>& v2)
{
    return Vertex<T>(v1.x + v2.x, v1.y + v2.y);
}

template <class T>
Vertex<T> operator*(const double k, const Vertex<T> &v1)
{
    return Vertex<T>(v1.x * k, v1.y * k);
}

#endif