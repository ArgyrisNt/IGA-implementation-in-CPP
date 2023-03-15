#ifndef H_BSPLINE_2D
#define H_BSPLINE_2D

#include "..\include\Bspline.h"
#include <iostream>
#include <vector>

class Bspline_2D
{
public:
    Bspline_2D() {}
    Bspline_2D(const Bspline &_bspline_x, const Bspline &_bspline_y)
        : bspline_x(_bspline_x), bspline_y(_bspline_y) {}

    Bspline &getBspline_x();
    Bspline &getBspline_y();
    KnotVector<double> &getKnotvector() { return bspline_x.getKnotvector(); }
    KnotVector<double> &YgetKnotvector() { return bspline_y.getKnotvector(); }

    void setBspline_x(Bspline &new_bspline_x);
    void setBspline_y(Bspline &new_bspline_y);

    int findSpanOfValue(double point);

private:
    Bspline bspline_x;
    Bspline bspline_y;
};

#include "..\src\Bspline_2D.cpp"

#endif