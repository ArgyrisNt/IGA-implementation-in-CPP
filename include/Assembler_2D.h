#pragma once

#include <iostream>
#include "..\include\Assembler.h"

class Assembler_2D : public Assembler
{
public:
    // Constructor
    Assembler_2D(double src, BoundCond &bcinfo, BsplineSurface &surface, std::vector<double> &trimm)
        : Assembler(src, bcinfo, surface.getBspline_x()), bspline_y(&surface.getBspline_y()),
          numberOfBasisFunctions(surface.getBspline_x().getNumberOfBasisFunctions() * surface.getBspline_y().getNumberOfBasisFunctions()), controlPoints(surface.getControlPoints()),
          trimming(trimm) {}

    Assembler_2D(double src, BoundCond &bcinfo, BsplineSurface &surface)
        : Assembler(src, bcinfo, surface.getBspline_x()), bspline_y(&surface.getBspline_y()),
          numberOfBasisFunctions(surface.getBspline_x().getNumberOfBasisFunctions() * surface.getBspline_y().getNumberOfBasisFunctions()), controlPoints(surface.getControlPoints()),
          trimming(std::vector<double>(3, 0.0)) {}

    // Destructor
    virtual ~Assembler_2D() {}

    // Member functions
    void assemble() override
    {
        computeTrimmedElements();
        computeStiffnessMatrix();
        computeRightHandSide();
        computeBoundary();
    }
    void plot_trimming();

    // Member getter functions
    Bspline& getBspline_y() { return *bspline_y; }
    std::vector<std::vector<double>>& getControlPoints() { return controlPoints; }
    const int getNumberOfBasisFunctions() const { return numberOfBasisFunctions; }

    // Member variables
    std::vector<double> trimming;
    std::vector<std::vector<std::pair<double, double>>> trimmed_triangles;

protected:
    // Functions for mapping coordinates
    std::vector<double> createTensorProduct(std::vector<double> &, std::vector<double> &);
    Matrix<double> Jacobian(double, double, std::pair<std::vector<double>, std::vector<double>> &, std::pair<std::vector<double>, std::vector<double>> &);
    
    std::pair<double, double> evaluate_trimming(double t);
    std::pair<double, double> evaluate_trimming_der(double t);
    double projection_on_trimming(std::pair<double, double> point);
    int checkIfOutside(std::pair<double, double> point);
    std::vector<std::vector<std::pair<double, double>>> construct_3_triangles(int id);
    std::vector<std::vector<std::pair<double, double>>> construct_2_triangles(int id);
    std::vector<std::vector<std::pair<double, double>>> construct_1_triangle(int id);
    std::pair<std::vector<double>, std::vector<double>> GaussPointsAndWeightsTria(double a, double b);
    int numberOfVoidNodesInElement(int elementX, int elementY);
    void categoriseElement(int elementX, int elementY);
    std::vector<std::vector<std::pair<double, double>>> divideElementInTriangles(int amount, int elementId);

    // Member local functions
    void computeTrimmedElements();
    void computeStiffnessMatrix();
    void computeRightHandSide();
    void computeBoundary();
    void computeTriangleStiffnessMatrix(std::vector<std::pair<double, double>> &triangle, int elementX, int elementY, Matrix<double> &A);
    void computeQuadStiffnessMatrix(int elementX, int elementY, Matrix<double> &A);
    void computeTriangleRightHandSide(std::vector<std::pair<double, double>> &triangle, int elementX, int elementY, std::vector<double> &b);
    void computeQuadRightHandSide(int elementX, int elementY, std::vector<double> &b);

    // Member variables
    Bspline *bspline_y;
    std::vector<std::vector<double>> controlPoints;
    const int numberOfBasisFunctions;
    std::vector<std::pair<bool, int>> trimmed_elements;
    std::vector<std::vector<std::pair<double, double>>> trimmed_elements_info; // contains all untrimmed vertices of each element
    std::vector<std::vector<std::pair<double, double>>> elements_vertices;     // contains all vertices of each element
};