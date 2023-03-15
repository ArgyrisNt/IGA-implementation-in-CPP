#ifndef H_ASSEMBLER_2D
#define H_ASSEMBLER_2D

#include <iostream>
#include "..\include\Assembler.h"
#include "..\include\Element.h"
#include "..\include\BsplineSurface.h"

class Assembler_2D : public Assembler<BsplineSurface>
{
public:
    Assembler_2D(double src, BoundCond &bcinfo, BsplineSurface &_surface)
        : Assembler<BsplineSurface>(src, bcinfo, _surface) {}

    virtual ~Assembler_2D() {}

    void assemble() override;

    int YspanOfValueInKnotVector(double value);
    
    void writeParameterSpaceToFile(std::string filename);
    void writeTrimmedTrianglesToFile(std::string filename);

    Bspline &getBspline_x();
    Bspline &getBspline_y();
    std::vector<std::vector<double>>& getControlPoints();
    const int getNumberOfBasisFunctions();
    double getDistinctKnotY(int position);
    std::vector<double> getDistinctKnotsY();
    TrimmingCurve &getTrimmingCurve();

    std::vector<Triangle<double>> trimmed_triangles;

protected:
    std::vector<int> computeActiveControlPoints(double g1, double g2);
    Matrix<double> Jacobian(std::vector<int> indices, std::vector<double> &dNxNy, std::vector<double> &NxdNy);

    void computeBoundary();
    int identifyBoundarySideOfBasisFunction(int i);
    void computeTrimmedElements();

    void computeStiffnessMatrixAndRighHandSide();
    void computeTriangleStiffnessMatrixAndRightHandSide(Triangle<double> &triangle, int elementX, int elementY);
    void computeQuadStiffnessMatrixAndRightHandSide(int elementX, int elementY);
    double computeStiffnessIntegral(int il_1, int jl_1, int il_2, int jl_2);
    double computeRightHandSideIntegral(int il_1, int il_2);

    std::vector<std::pair<double, double>> YGaussPointsAndWeights;
    std::vector<Element> elements;
};

#include "..\src\Assembler_2D.cpp"

#endif