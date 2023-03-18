#include "..\include\MultiBspline.h"
#include <iostream>
#include <string>


Bspline &MultiBspline::getBspline(int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    return bsplines[dim];
}

int MultiBspline::getDegree(int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    return bsplines[dim].getDegree();
}

KnotVector<double> &MultiBspline::getKnotvector(int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    return bsplines[dim].getKnotvector();
}

int MultiBspline::getDimension()
{
    return dimension;
}



void MultiBspline::setBspline(const Bspline &new_bspline, int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    bsplines[dim] = new_bspline;
}



int MultiBspline::findSpanOfValue(double point, int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    return bsplines[dim].findSpanOfValue(point);
}

void MultiBspline::plot2D(int resolution, const std::string &filename)
{
    for (int i = 0; i < dimension; ++i)
        bsplines[i].plot2D(resolution, std::to_string(i + 1) + filename);
}