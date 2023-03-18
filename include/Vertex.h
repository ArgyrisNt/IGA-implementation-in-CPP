#include <iostream>

template <class T>
class Vertex
{
public:
    Vertex() {}
    Vertex(T _x, T _y) : x(_x), y(_y) {}

    void set(T _x, T _y) { x = _x; y = _y; }

    bool operator==(Vertex<T>& vertex2)
    {
        return (fabs(x - vertex2.x) < 1e-7 && fabs(y - vertex2.y) < 1e-7);
    }

    bool operator!=(Vertex<T>& vertex2)
    {
        return !((*this) == vertex2);
    }


    T x;
    T y;
};