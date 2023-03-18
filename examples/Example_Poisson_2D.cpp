// Solve Poisson 2D problem on a Bspline quarter annulus

#include <iostream>
#include "..\include\Assembler_2D.h"
#include "..\include\Poisson.h"
#include "..\include\BsplineSurface.h"

int main()
{   
    // - - - - - B-spline basis on x-direction - - - - -
    int x_degree = 2;
    std::vector<double> x_values{ 0.0,0.0,0.0,1.0,1.0,1.0 };
    std::vector<double> x_weights{1.0, 1.0, 1.0};
    KnotVector<double> x_knotVector(x_degree, x_values, x_weights);
    Bspline bspline_x(x_knotVector);

    // - - - - - B-spline basis on y-direction - - - - - 
    int y_degree = 2;
    std::vector<double> y_values{ 0.0,0.0,0.0,0.5,1.0,1.0,1.0 };
    std::vector<double> y_weights{ 1.0,1.0,1.0,1.0 };
    KnotVector<double> y_knotVector(y_degree, y_values, y_weights);
    Bspline bspline_y(y_knotVector);

    // - - - - - B-spline surface - - - - -
    std::vector<std::vector<double>> controlPoints{ {0.0, 0.0}, {0.0, 4.0}, {4.0, 8.0}, {8.0, 8.0},
                                                    {2.0, 0.0}, {2.0, 3.0}, {5.0, 6.0}, {8.0, 6.0},
                                                    {4.0, 0.0}, {4.0, 2.0}, {6.0, 4.0}, {8.0, 4.0} };
    TrimmingCurve trimmingCurve(Vertex<double>(0.0, 0.0), 0.0);
    BsplineSurface surface(std::vector<Bspline>{bspline_x, bspline_y}, controlPoints, trimmingCurve);
    for (int i = 0; i < 1; ++i)
    {
        surface.uniformRefine_x();
        surface.uniformRefine_y();
    }

    // - - - - - Assempler info - - - - -
    double src = 3.0;
    std::pair<std::string, double> west = std::make_pair("Dirichlet", 0.0);
    std::pair<std::string, double> east = std::make_pair("Dirichlet", 0.0);
    std::pair<std::string, double> north = std::make_pair("Dirichlet", 0.0);
    std::pair<std::string, double> south = std::make_pair("Dirichlet", 0.0);
    BoundCond boundaryConditions(west, east, north, south);
    Assembler_2D assembler(src, boundaryConditions, surface);
    assembler.assemble();

    // - - - - - Enforce boundary conditions - - - - -
    std::string mode("Ellimination");
    assembler.enforceBoundaryConditions(mode);

    // - - - - - Poisson info - - - - -
    Poisson<Assembler_2D> poisson(assembler, Solver::SOR);
    int iterations = 50;
    poisson.solve(iterations);
    std::cout << poisson.getSolution();

    // - - - - - Write solution data - - - - - 
    int resolution = 100;
    //surface.plot2D(resolution, "surface.dat");
    surface.plot3D(resolution, poisson.getSolution(), "solution.dat");
    writeParameterSpaceToFile(assembler, "parameterSpace.dat");

    return 0;
}