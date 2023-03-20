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
    std::vector<Vertex<double>> controlPoints{ {0.0, 0.0}, {0.0, 3.0}, {0.0, 6.0},
                                              {3.0, 0.0}, {3.0, 3.0}, {3.0, 6.0},
                                              {6.0, 0.0}, {6.0, 3.0}, {6.0, 6.0} };
    TrimmingCurve trimmingCurve(Vertex<double>(1.0, 0.0), 0.2 /*std::make_pair(1.8,0.5), 1.0*/);
    trimmingCurve.plot();
    BsplineSurface surface(std::vector<Bspline>{bspline_x, bspline_y}, controlPoints, trimmingCurve);
    for (int i = 0; i < 2; ++i) surface.uniformRefine_x();
    for (int i = 0; i < 2; ++i) surface.uniformRefine_y();

    // - - - - - Assempler info - - - - -
    double src = 3.0;
    std::pair<std::string, double> west = std::make_pair("Dirichlet", 0.0);
    std::pair<std::string, double> east = std::make_pair("Neumann", 0.0);
    std::pair<std::string, double> north = std::make_pair("Dirichlet", 0.0);
    std::pair<std::string, double> south = std::make_pair("Neumann", 0.0);
    BoundCond boundaryConditions(west, east, north, south);
    Assembler_2D assembler(src, boundaryConditions, surface);
    assembler.assemble();
    std::string mode("Ellimination");
    assembler.enforceBoundaryConditions(mode);

    // - - - - - Poisson info - - - - -
    Poisson<Assembler_2D> poisson(assembler, Solver::GaussSeidel);
    poisson.solve();
    std::cout << poisson.getSolution();

    // - - - - - Plot solution - - - - - 
    int resolution = 100;
    surface.plotVectorOnEntity(resolution, poisson.getSolution(), "solution.dat");
    writeParameterSpaceToFile(assembler, "parameterSpace.dat");
    writeTrimmedTrianglesToFile(assembler, "trimmedTriangles.obj");

    return 0;
}