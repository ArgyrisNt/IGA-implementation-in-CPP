#ifndef H_MULTIBSPLINE
#define H_MULTIBSPLINE

#include "Bspline.h"
#include <iostream>
#include <vector>

class MultiBspline
{
public:
    MultiBspline() {}
    MultiBspline(const std::vector<Bspline> &newBsplines) : bsplines(newBsplines), dimension(bsplines.size()) {}

    ~MultiBspline() {}

    Bspline &getBspline(int dim);
    KnotVector<double> &getKnotvector(int dim);
    const int getDegree(int dim) const;
    const int getDimension() const 
    { return dimension; }

    int findSpanOfValue(double point, int dim);
    void plot(int resolution, const std::string &filename);
    void setBspline(const Bspline &new_bspline, int dim);

private:
    std::vector<Bspline> bsplines;
    int dimension;
};

inline Bspline &MultiBspline::getBspline(int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    return bsplines[dim];
}

inline const int MultiBspline::getDegree(int dim) const
{
    assert(dim >= 0 && dim < bsplines.size());
    return bsplines[dim].getDegree();
}

inline KnotVector<double> &MultiBspline::getKnotvector(int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    return bsplines[dim].getKnotvector();
}

inline void MultiBspline::setBspline(const Bspline &new_bspline, int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    bsplines[dim] = new_bspline;
}

inline int MultiBspline::findSpanOfValue(double point, int dim)
{
    assert(dim >= 0 && dim < bsplines.size());
    return bsplines[dim].findSpanOfValue(point);
}

#include "..\src\MultiBspline.cpp"

#endif