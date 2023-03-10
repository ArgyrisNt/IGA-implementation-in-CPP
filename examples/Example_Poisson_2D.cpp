// Solve Poisson 2D problem on a Bspline quarter annulus

#include <iostream>
#include "..\IGA.h"

int main()
{   
    // - - - - - B-spline basis on x-direction - - - - -
    int x_degree = 2;
    std::vector<double> x_values{ 0.0,0.0,0.0,1.0,1.0,1.0 };
    KnotVector<double> x_knotVector(x_degree, x_values);
    std::vector<double> x_weights{ 1.0,1.0,1.0 };
    Bspline bspline_x(x_degree, x_knotVector, x_weights);

    // - - - - - B-spline basis on y-direction - - - - - 
    int y_degree = 2;
    std::vector<double> y_values{ 0.0,0.0,0.0,0.5,1.0,1.0,1.0 };
    std::vector<double> y_weights{ 1.0,1.0,1.0,1.0 };
    KnotVector<double> y_knotVector(y_degree, y_values);
    Bspline bspline_y(y_degree, y_knotVector, y_weights);

    // - - - - - B-spline surface - - - - -
    std::vector<std::vector<double>> controlPoints{ {0.0, 0.0}, {0.0, 4.0}, {4.0, 8.0}, {8.0, 8.0},
                                                    {2.0, 0.0}, {2.0, 3.0}, {5.0, 6.0}, {8.0, 6.0},
                                                    {4.0, 0.0}, {4.0, 2.0}, {6.0, 4.0}, {8.0, 4.0} };
    // std::vector<std::vector<double>> controlPoints{ {0.0, 0.0}, {0.0, 2.0}, {0.0, 4.0},
    //                                                 {2.0, 0.0}, {2.0, 2.0}, {2.0, 4.0},
    //                                                 {4.0, 0.0}, {4.0, 2.0}, {4.0, 4.0} };
    BsplineSurface surface(bspline_x, bspline_y, controlPoints);
    for (int i = 0; i < 1; i++)
    {
        surface.uniformRefine_x();
        surface.uniformRefine_y();
    }

    // - - - - - Assempler info - - - - -
    double src = 3.0;
    BoundCond _bc("Dirichlet", "Dirichlet", "Dirichlet", "Dirichlet", 0.0, 0.0, 0.0, 0.0); // left-right-top-bottom
    TrimmingCurve trimmingCurve(std::make_pair(0.0, 0.0), 0.0);
    Assembler_2D ass2(src, _bc, surface, trimmingCurve);
    ass2.assemble();

    // - - - - - Enforce boundary conditions - - - - -
    std::string mode("Ellimination");
    ass2.enforceBoundaryConditions(mode);

    // - - - - - Poisson info - - - - -
    Poisson_2D poisson(ass2, Solver::SOR);
    poisson.setSolution(poisson.getSolver()->solve(ass2.getStiffnessMatrix(), ass2.getRightHandSide(), 50));
    poisson.expandSolutionOnBoundary();
    std::cout << poisson.getSolution();

    // - - - - - Write solution data - - - - - 
    poisson.plotSolution(100); // resolution = 100

    return 0;
}