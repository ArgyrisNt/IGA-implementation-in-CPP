// Solve Diffusion 1D problem on a Bspline line

#include <iostream>
#include <string>
#include "..\IGA.h"

double init_cond(double val)
{
    double res;
    if (val >= 0 && val <= 0.5)
    {
        return 2 * val;
    }
    else if (val > 0.5 && val <= 1.0)
    {
        return 2 - 2 * val;
    }
    // if (val >=0 && val <= 1)
    // {
    // 	return exp(-pow((val - 0.5), 2) / (pow(0.05, 2)));
    // }
    else
    {
        std::cout << "Invalid input for function." << std::endl;
        throw std::invalid_argument("Invalid input");
    }
}

int main()
{   
    // - - - - - B-spline basis on x-direction - - - - - 
    double start = 0.0;
    double end = 1.0;
    int p = 3;
    int numElems = 200;
    std::vector<double> W{};
    std::vector<std::vector<double>> ctrlPts;
    Bspline bspline_x(start, end, p, numElems, W);
    for (int i = 0; i < numElems + p; i++)
    {
        ctrlPts.push_back({(1.0 / (numElems + p - 1)) * (double)(i), 0.0});
        W.push_back(1.0);
    }
    bspline_x.setWeights(W);

    // - - - - - B-spline curve - - - - -
    BsplineCurve curve(bspline_x, ctrlPts);

    // - - - - - Assempler info - - - - -
    double src = 0.0;
    BoundCond _bc("Dirichlet", "Dirichlet", 0.0, 0.0);
    double coef = 1.0; // coefficient of diffusion problem
    double t_end = 0.05; // 1
    double numSteps = 10.0; // 200
    double Dt = t_end / numSteps; // time step
    DiffusionAssembler_1D ass(src, _bc, curve, coef, Dt);
    ass.assemble();

    // - - - - - Enforce boundary conditions - - - - -
    std::string mode("Ellimination");
    ass.enforceBoundary(mode);

    // - - - - - Apply and plot initial condition - - - - -
    Poisson<DiffusionAssembler_1D> diffusion(ass, Solver::GaussSeidel);
    std::vector<double> init_sol = ass.applyInitCond(init_cond);
    diffusion.setSolution(init_sol);
    diffusion.plotSol("curve.dat", "0sol.dat");

    // - - - - - Solve - - - - - 
    for (int t = 0; t < numSteps; t++) // until 6
    {
        std::cout << std::endl << "---------------- " << t + 1 << " step ----------------";
        std::vector<double> b = ass.nextStep(diffusion.getSolution()); // build next rhs
        diffusion.setSolution(diffusion.getSolver()->solve(ass.getSysMat(), b));
        diffusion.constructSol();
        diffusion.plotSol("curve.dat", std::to_string(t + 1) + "sol.dat");
    }

    return 0;
}