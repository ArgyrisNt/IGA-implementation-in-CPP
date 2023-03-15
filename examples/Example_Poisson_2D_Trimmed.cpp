// Solve Poisson 2D problem on a Bspline trimmed square domain

#include <iostream>
#include "..\include\Assembler_2D.h"
#include "..\include\Poisson.h"

int main()
{   
    // - - - - - B-spline basis on x-direction - - - - - 
    int x_degree = 2;
    std::vector<double> x_values{ 0.0,0.0,0.0,1.0,1.0,1.0 };
    std::vector<double> x_weights{1.0, 1.0, 1.0};
    KnotVector<double> knotVector(x_degree, x_values, x_weights);
    Bspline bspline_x(knotVector);

    // - - - - - B-spline basis on y-direction - - - - - 
    Bspline bspline_y(bspline_x);

    // - - - - - B-spline surface - - - - -
    std::vector<std::vector<double>> controlPoints{ {0.0, 0.0}, {0.0, 3.0}, {0.0, 6.0},
                                              {3.0, 0.0}, {3.0, 3.0}, {3.0, 6.0},
                                              {6.0, 0.0}, {6.0, 3.0}, {6.0, 6.0} };
    TrimmingCurve trimmingCurve(Vertex<double>(1.0, 0.0), 0.2 /*std::make_pair(1.8,0.5), 1.0*/);
    trimmingCurve.plot();
    BsplineSurface surface(bspline_x, bspline_y, controlPoints, trimmingCurve);
    for (int i = 0; i < 2; i++) surface.uniformRefine_x();
    for (int i = 0; i < 2; i++) surface.uniformRefine_y();

    // - - - - - Assempler info - - - - -
    double src = 3.0;
    BoundCond _bc("Dirichlet", "Neumann", "Dirichlet", "Neumann", 0.0, 0.0, 0.0, 0.0); // left-right-top-bottom
    Assembler_2D ass2(src, _bc, surface);
    ass2.assemble();

    // - - - - - Enforce boundary conditions - - - - -
    std::string mode("Ellimination");
    ass2.enforceBoundaryConditions(mode);

    // - - - - - Poisson info - - - - -
    Poisson<Assembler_2D> poisson(ass2, Solver::GaussSeidel);
    poisson.setSolution(poisson.getSolver()->solve());
    std::cout << poisson.getSolution();

    // - - - - - Write solution data - - - - - 
    int resolution = 250;
    surface.plot3D(resolution, poisson.getSolution(), "solution.dat");
    ass2.writeParameterSpaceToFile("parameterSpace.dat");
    ass2.writeTrimmedTrianglesToFile("trimmedTriangles.obj");

    return 0;
}