// Solve Poisson 2D problem on a Bspline trimmed square domain

#include <iostream>
#include "..\IGA.h"

int main()
{   
    // - - - - - B-spline basis on x-direction - - - - - 
    int x_degree = 2;
    std::vector<double> x_values{ 0.0,0.0,0.0,1.0,1.0,1.0 };
    KnotVector<double> knotVector(x_degree, x_values);
    std::vector<double> x_weights{ 1.0,1.0,1.0 };
    Bspline bspline_x(x_degree, knotVector, x_weights);

    // - - - - - B-spline basis on y-direction - - - - - 
    Bspline bspline_y(bspline_x);

    // - - - - - B-spline surface - - - - -
    std::vector<std::vector<double>> controlPoints{ {0.0, 0.0}, {0.0, 3.0}, {0.0, 6.0},
                                              {3.0, 0.0}, {3.0, 3.0}, {3.0, 6.0},
                                              {6.0, 0.0}, {6.0, 3.0}, {6.0, 6.0} };
    BsplineSurface surface(bspline_x, bspline_y, controlPoints);
    for (int i = 0; i < 1; i++) surface.uniformRefine_x();
    for (int i = 0; i < 1; i++) surface.uniformRefine_y();

    // - - - - - Assempler info - - - - -
    double src = 3.0;
    BoundCond _bc("Dirichlet", "Neumann", "Dirichlet", "Dirichlet", 0.0, 0.0, 0.0, 0.0); // left-right-top-bottom
    std::vector<double> trimming_curve{1.7, 0.5, 1.0}; // circle with center = (1.7,0.5) and radius = 1.0;
    Assembler_2D ass2(src, _bc, surface, trimming_curve);
    ass2.plot_trimming();
    ass2.assemble();

    // - - - - - Enforce boundary conditions - - - - -
    std::string mode("Ellimination");
    ass2.enforceBoundaryConditions(mode);

    // - - - - - Poisson info - - - - -
    Poisson_2D poisson(ass2, Solver::GaussSeidel);
    poisson.setSolution(poisson.getSolver()->solve(ass2.getStiffnessMatrix(), ass2.getRightHandSide()));
    poisson.expandSolutionOnBoundary();
    std::cout << poisson.getSolution();

    // - - - - - Write solution data - - - - - 
    poisson.plotSolution(100); // resolution = 500

    return 0;
}