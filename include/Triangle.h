#include <iostream>
#include "Vertex.h"

template <class T>
class Triangle
{
public:
    Triangle(Vertex<T> &_vertex1, Vertex<T> &_vertex2, Vertex<T> &_vertex3)
        : vertex1(_vertex1), vertex2(_vertex2), vertex3(_vertex3) {}

    ~Triangle() {}

    Vertex<T> vertex1;
    Vertex<T> vertex2;
    Vertex<T> vertex3;
};