#ifndef H_MULTIBSPLINE
#define H_MULTIBSPLINE

#include "..\include\Bspline.h"
#include <iostream>
#include <vector>

class MultiBspline
{
public:
    MultiBspline() {}
    MultiBspline(const std::vector<Bspline> &newBsplines) : bsplines(newBsplines), dimension(bsplines.size()) {}

    virtual ~MultiBspline() {}

    Bspline &getBspline(int dim);
    KnotVector<double> &getKnotvector(int dim);
    int getDegree(int dim);
    int getDimension();

    void setBspline(const Bspline &new_bspline, int dim);

    int findSpanOfValue(double point, int dim);
    virtual void plot2D(int resolution, const std::string &filename);

private:
    std::vector<Bspline> bsplines;
    int dimension;
};

#include "..\src\MultiBspline.cpp"

#endif