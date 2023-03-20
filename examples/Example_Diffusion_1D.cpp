// Solve Diffusion 1D problem on a Bspline line

#include <iostream>
#include <string>
#include "..\include\DiffusionAssembler_1D.h"
#include "..\include\Poisson.h"

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
    // if (val >=0.0 && val <= 1.0)
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
    // - - - - - B-spline basis - - - - - 
    double start = 0.0;
    double end = 1.0;
    int degree = 3;
    int numberOfElements = 200;
    std::vector<double> weights{};
    KnotVector<double> knotVector(start, end, degree, numberOfElements, weights);
    std::vector<Vertex<double>> controlPoints;
    for (int i = 0; i < numberOfElements + degree; ++i)
    {
        controlPoints.push_back({(1.0 / (numberOfElements + degree - 1)) * (double)(i), 0.0});
        weights.push_back(1.0);
    }
    knotVector.setWeights(weights);
    Bspline bspline_x(knotVector);

    // - - - - - B-spline curve - - - - -
    BsplineCurve curve(std::vector<Bspline>{bspline_x}, controlPoints);

    // - - - - - Assempler info - - - - -
    double src = 0.0;
    std::pair<std::string, double> west = std::make_pair("Dirichlet", 0.0);
    std::pair<std::string, double> east = std::make_pair("Dirichlet", 0.0);
    BoundCond boundaryConditions(west, east);
    double coef = 1.0; // coefficient of diffusion problem
    double t_end = 0.05;
    double numSteps = 10.0;
    double Dt = t_end / numSteps; // time step
    DiffusionAssembler_1D ass(src, boundaryConditions, curve, coef, Dt);
    ass.assemble();
    std::string mode("Ellimination");
    ass.enforceBoundaryConditions(mode);

    // - - - - - Apply and plot initial condition - - - - -
    Poisson<DiffusionAssembler_1D> diffusion(ass, Solver::GaussSeidel);
    std::vector<double> init_sol = ass.applyInitialCondition(init_cond);
    diffusion.applyInitialCondition(init_sol);
    int resolution = 100;
    curve.plotVectorOnEntity(resolution, diffusion.getSolution(), "0solution.dat");

    // - - - - - Solve and plot solution on every step - - - - - 
    for (int t = 0; t < numSteps; ++t)
    {
        std::cout << std::endl << "---------------- " << t + 1 << " step ----------------";
        std::vector<double> b = ass.nextStep(diffusion.getSolution()); // build next rhs
        diffusion.updateRhs(b);
        diffusion.solve(20);
        curve.plotVectorOnEntity(resolution, diffusion.getSolution(), std::to_string(t + 1) + "solution.dat");
    }

    return 0;
}