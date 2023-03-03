#pragma once

#include <iostream>
#include "..\include\Matrix.h"
#include "..\include\BoundCond.h"
#include "..\include\Bspline.h"

class Assembler
{
public:
    // Cunstructor
    Assembler(double src, BoundCond& bcinfo, Bspline& bas) : f(src), bc(&bcinfo), bspline_x(&bas) {}

    // Destructor
    virtual ~Assembler() {}

    // Member functions
    virtual void assemble() = 0;
    virtual void applyBoundEllimination();
    virtual void applyBoundMultipliers();
    virtual void enforceBoundary(std::string&);

    // Member getter functions
    Matrix<double>& getStiff() { return stiff; }
    std::vector<double>& getRhs() { return rhs; }
    Bspline& getBspline_x() { return *bspline_x; }
    BoundCond& getBc() { return *bc; }  
    std::string& getBoundaryMode() { return boundaryMode; }
    std::vector<std::pair<int, int>>& getBoundaryIds() { return boundary_ids; }
    
protected:
    // Member variables
    Matrix<double> stiff;
    std::vector<double> rhs;
    Bspline* bspline_x;
    double f;
    BoundCond* bc;
    std::string boundaryMode;
    std::vector<std::pair<int, int>> boundary_ids;
};