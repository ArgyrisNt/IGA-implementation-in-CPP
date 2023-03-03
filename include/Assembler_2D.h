#pragma once

#include <iostream>
#include "..\include\Assembler.h"

class Assembler_2D : public Assembler
{
public:
    // Constructor
    Assembler_2D(double src, BoundCond &bcinfo, BsplineSurface &surface, std::vector<double>& trimm)
        : Assembler(src, bcinfo, surface.getBspline_x()), bspline_y(&surface.getBspline_y()),
            nOF(surface.getBspline_x().getNOF() * surface.getBspline_y().getNOF()), ctrlPts(surface.getCtrlPts()),
            trimming(trimm) {}

    Assembler_2D(double src, BoundCond &bcinfo, BsplineSurface &surface)
        : Assembler(src, bcinfo, surface.getBspline_x()), bspline_y(&surface.getBspline_y()),
          nOF(surface.getBspline_x().getNOF() * surface.getBspline_y().getNOF()), ctrlPts(surface.getCtrlPts()),
          trimming(std::vector<double>(3,0.0)) {}

    // Destructor
    virtual ~Assembler_2D() {}

    // Member functions
    void assemble() override
    {
        calcTrimmed();
        calcStiff();
        calcRhs();
        calcBound();
    }

    std::vector<double> createTensorProduct(std::vector<double>&, std::vector<double>&);
    Matrix<double> calcJacobian(double, double, int, int, std::vector<double>&, std::vector<double>&);
    std::pair<std::vector<double>,std::vector<double>> Map2Physical(Matrix<double>&, std::vector<double>&, std::vector<double>&);

    // Member getter functions
    Bspline& getBspline_y() { return *bspline_y; }
    std::vector<std::vector<double>>& getCtrlPts() { return ctrlPts; }
    const int getNOF() const { return nOF; } 

    void plot_trimming();
    std::vector<double> trimming;
    std::vector<std::vector<std::pair<double, double>>> trimmed_triangles;

protected:
    std::vector<std::pair<bool, int>> trimmed_elements;
    std::vector<std::vector<std::pair<double, double>>> trimmed_elements_info; // contains all untrimmed vertices of each element
    std::vector<std::vector<std::pair<double, double>>> elements_vertices; // contains all vertices of each element
    std::pair<double, double> evaluate_trimming(double t);
    std::pair<double, double> evaluate_trimming_der(double t);
    double projection_on_trimming(std::pair<double, double> point);
    int checkIfOutside(std::pair<double, double> point);
    std::vector<std::vector<std::pair<double, double>>> construct_3_triangles(int id);
    std::vector<std::vector<std::pair<double, double>>> construct_2_triangles(int id);
    std::vector<std::vector<std::pair<double, double>>> construct_1_triangle(int id);
    std::pair<std::vector<double>, std::vector<double>> calcGaussPtsTria(double a, double b);

    // Member local functions
    void calcTrimmed();
    void calcStiff();
    void calcRhs();
    void calcBound();

    // Member variables
    Bspline* bspline_y;
    std::vector<std::vector<double>> ctrlPts;
    const int nOF;
};