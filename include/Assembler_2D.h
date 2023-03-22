#ifndef H_ASSEMBLER_2D
#define H_ASSEMBLER_2D

#include <iostream>
#include "Assembler.h"
#include "Element.h"
#include "BsplineSurface.h"

class Assembler_2D : public Assembler<BsplineSurface>
{
public:
    Assembler_2D(const double src, const BoundCond &bcinfo, const BsplineSurface &_surface)
        : Assembler<BsplineSurface>(src, bcinfo, _surface) {}

    ~Assembler_2D() {}

    Bspline &getBspline_x() 
    { return bsplineEntity->getMultiBspline().getBspline(0); }

    Bspline &getBspline_y() 
    { return bsplineEntity->getMultiBspline().getBspline(1); }

    const int getNumberOfBasisFunctions()
    { return getBspline_x().getNumberOfBasisFunctions() * getBspline_y().getNumberOfBasisFunctions(); }

    const TrimmingCurve &getTrimmingCurve() const 
    { return bsplineEntity->getTrimmingCurve(); }

    const std::vector<Triangle<double>> &getTrimmedTriangles() const 
    { return trimmed_triangles; }

private:
    std::vector<int> computeActiveControlPoints(double g1, double g2);
    Matrix<double> Jacobian(const std::vector<int> &indices, const std::vector<double> &dNxNy, const std::vector<double> &NxdNy);

    void computeBoundary() override;
    int identifyBoundarySideOfBasisFunction(int i);
    void computeTrimmedElements();

    void computeStiffnessMatrixAndRightHandSide() override;
    void AssistComputeStiffnessMatrixAndRightHandSide(int spanX, int spanY);
    double computeStiffnessIntegral(int il_1, int jl_1, int il_2, int jl_2);
    double computeRightHandSideIntegral(int il_1, int il_2);

    std::vector<std::pair<double, double>> YGaussPointsAndWeights;
    std::vector<Element> elements;
    std::vector<Triangle<double>> trimmed_triangles;
};

void writeParameterSpaceToFile(Assembler_2D& assembler, const std::string &filename)
{
    std::ofstream my_file(filename);
	my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
	my_file << "zone t= " << "\"1\"" << ",i=" << assembler.getDistinctKnots(0).size() << ",j=" << assembler.getDistinctKnots(1).size() << "\n";
    for (int j = 0; j < assembler.getDistinctKnots(1).size(); ++j)
    {
        for (int i = 0; i < assembler.getDistinctKnots(0).size(); ++i)
        {
            my_file << assembler.getDistinctKnots(0)[i] << " " << assembler.getDistinctKnots(1)[j] << "\n";
        }
	}
	my_file.close();
}

void writeTrimmedTrianglesToFile(Assembler_2D &assembler, const std::string &filename)
{
    std::ofstream my_file(filename);
    my_file << "# " << assembler.getTrimmedTriangles().size() << "\n";
    int cnt = 0;
    for (auto triangle : assembler.getTrimmedTriangles())
    {
        cnt += 2;
        my_file << "v " << triangle.vertex1.x << " " << triangle.vertex1.y << "\n";
        my_file << "v " << triangle.vertex2.x << " " << triangle.vertex2.y << "\n";
        my_file << "l " << cnt - 1 << " " << cnt << "\n";
        cnt += 2;
        my_file << "v " << triangle.vertex2.x << " " << triangle.vertex2.y << "\n";
        my_file << "v " << triangle.vertex3.x << " " << triangle.vertex3.y << "\n";
        my_file << "l " << cnt - 1 << " " << cnt << "\n";
        cnt += 2;
        my_file << "v " << triangle.vertex3.x << " " << triangle.vertex3.y << "\n";
        my_file << "v " << triangle.vertex1.x << " " << triangle.vertex1.y << "\n";
        my_file << "l " << cnt - 1 << " " << cnt << "\n";
    }
    my_file.close();
}

#include "..\src\Assembler_2D.cpp"

#endif