#include "..\include\Bspline_2D.h"
#include <iostream>

Bspline &Bspline_2D::getBspline_x()
{
    return bspline_x;
}

Bspline &Bspline_2D::getBspline_y()
{
    return bspline_y;
}

void Bspline_2D::setBspline_x(Bspline &new_bspline_x)
{
    bspline_x = new_bspline_x;
}

void Bspline_2D::setBspline_y(Bspline &new_bspline_y)
{
    bspline_y = new_bspline_y;
}

int Bspline_2D::findSpanOfValue(double value)
{
    return bspline_x.findSpanOfValue(value);
}
