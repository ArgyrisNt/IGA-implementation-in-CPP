// Solve Poisson 2D problem on a Bspline quarter annulus

#include <iostream>
#include "..\IGA.h"

int main()
{   
    // - - - - - B-spline basis on x-direction - - - - - 
    int p = 2;
    std::vector<double> U{ 0.0,0.0,0.0,1.0,1.0,1.0 };
    std::vector<double> W{ 1.0,1.0,1.0 };
    Bspline bspline_x(p, U, W);

    // - - - - - B-spline basis on y-direction - - - - - 
    int q = 2;
    std::vector<double> V{ 0.0,0.0,0.0,0.5,1.0,1.0,1.0 };
    std::vector<double> Q{ 1.0,1.0,1.0,1.0 };
    Bspline bspline_y(q, V, Q);

    // - - - - - B-spline surface - - - - -
    std::vector<std::vector<double>> ctrlPts{ {0.0, 0.0}, {0.0, 4.0}, {4.0, 8.0}, {8.0, 8.0},
                                              {2.0, 0.0}, {2.0, 3.0}, {5.0, 6.0}, {8.0, 6.0},
                                              {4.0, 0.0}, {4.0, 2.0}, {6.0, 4.0}, {8.0, 4.0} };
    BsplineSurface surface(bspline_x, bspline_y, ctrlPts);
    for (int i = 0; i < 0; i++)
    {
        surface.uniformRefine_y();
    }

    // - - - - - Assempler info - - - - -
    double src = 3.0;
    BoundCond _bc("Dirichlet", "Dirichlet", "Neumann", "Dirichlet", 0.0, 0.0, 0.0, 0.0); // left-right-top-bottom
    Assembler_2D ass2(src, _bc, surface);
    ass2.assemble();

    // - - - - - Enforce boundary conditions - - - - -
    std::string mode("Ellimination");
    ass2.enforceBoundary(mode);

    // - - - - - Poisson info - - - - -
    Poisson_2D poisson(ass2, Solver::ConjugateGradient);
    poisson.setSolution(poisson.getSolver()->solve(ass2.getStiff(), ass2.getRhs(), 50));
    poisson.constructSol();
    std::cout << poisson.getSolution();

    // - - - - - Write solution data - - - - - 
    poisson.plotSol(100); // resolution = 100

    return 0;
}