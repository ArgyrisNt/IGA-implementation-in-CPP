#pragma once

#include <iostream>
#include "..\include\Assembler.h"
#include "..\include\Element.h"

class Assembler_2D : public Assembler
{
public:
    Assembler_2D(double src, BoundCond &bcinfo, BsplineSurface &surface, TrimmingCurve& _trimmingCurve)
        : Assembler(src, bcinfo, surface.getBspline_x()), bspline_y(&surface.getBspline_y()),
          numberOfBasisFunctions(surface.getBspline_x().getNumberOfBasisFunctions() * surface.getBspline_y().getNumberOfBasisFunctions()),
        controlPoints(surface.getControlPoints()), trimmingCurve(_trimmingCurve) {}

    virtual ~Assembler_2D() {}

    void assemble() override;

    int YspanOfValueInKnotVector(double value);
    
    void writeParameterSpaceToFile(std::string filename);

    Bspline& getBspline_y() { return *bspline_y; }
    std::vector<std::vector<double>>& getControlPoints() { return controlPoints; }
    const int getNumberOfBasisFunctions() const { return numberOfBasisFunctions; }

    std::vector<Triangle<double>> trimmed_triangles;
    TrimmingCurve trimmingCurve;

protected:
    std::vector<double> createTensorProduct(std::vector<double> &, std::vector<double> &);
    std::vector<int> computeActiveControlPoints(double g1, double g2);
    Matrix<double> Jacobian(double, double, std::pair<std::vector<double>, std::vector<double>> &, std::pair<std::vector<double>, std::vector<double>> &);

    std::vector<std::pair<double, double>> GaussPointsAndWeightsTria(double a, double b);

    bool basisFunctionHasNotAlreadyMarked(int id);

    void YcomputeDistinctKnots();
    void computeBoundary();
    void computeTrimmedElements();

    void computeStiffnessMatrixAndRighHandSide();
    void computeTriangleStiffnessMatrixAndRightHandSide(Triangle<double> &triangle, int elementX, int elementY, Matrix<double> &A, std::vector<double> &b);
    void computeQuadStiffnessMatrixAndRightHandSide(int elementX, int elementY, Matrix<double> &A, std::vector<double>& b);
    double computeStiffnessIntegral(std::vector<std::pair<double, double>> &gauss_x, int il_1, int jl_1, std::vector<std::pair<double, double>> &gauss_y, int il_2, int jl_2);
    double computeRightHandSideIntegral(std::vector<std::pair<double,double>> &gauss_x, int il_1, std::vector<std::pair<double, double>> &gauss_y, int il_2);

    std::vector<double> YdistinctKnots;
    Bspline *bspline_y;
    std::vector<std::vector<double>> controlPoints;
    const int numberOfBasisFunctions;
    std::vector<Element> elements;
};