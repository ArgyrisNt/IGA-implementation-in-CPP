#include <iostream>
#include "..\include\Poisson.h"

template <class T>
Poisson<T>::Poisson(T &newAssembler, Solver &newSolver)
{
    assembler = &newAssembler;
    solver = &newSolver;
    solver->setLeftAndRightHandSides(assembler->getSystemMatrix(), assembler->getRightHandSide());
}

template <class T>
void Poisson<T>::expandSolutionOnBoundary()
{
    std::vector<double> newSolution;
    if (assembler->boundaryMode == "Multipliers")
    {
        solution.erase(solution.begin(), solution.begin() + assembler->boundaryBasisFunctions.size());
    }
    else
    {
	    int j = 0;
	    for (int i = 0; i < assembler->getNumberOfBasisFunctions(); i++)
	    {
            auto it = std::find_if(assembler->boundaryBasisFunctions.begin(), assembler->boundaryBasisFunctions.end(), CompareFirst(i));
            if (it != assembler->boundaryBasisFunctions.end())
            {
                int position = it - assembler->boundaryBasisFunctions.begin();
                if (assembler->boundaryBasisFunctions[position].second == 1)
                {
                    newSolution.push_back(assembler->boundaryConditions->getWestValue());
                }
                else if (assembler->boundaryBasisFunctions[position].second == 2)
                {
                    newSolution.push_back(assembler->boundaryConditions->getEastValue());
                }
                else if (assembler->boundaryBasisFunctions[position].second == 3)
                {
                    newSolution.push_back(assembler->boundaryConditions->getSouthValue());
                }
                else if (assembler->boundaryBasisFunctions[position].second == 4)
                {
                    newSolution.push_back(assembler->boundaryConditions->getNorthValue());
                }							
		    }
		    else
		    {
                newSolution.push_back(solution[j]);
                j++;
		    }
	    }
        solution = newSolution;
    }
}

template <class T>
std::vector<double> Poisson<T>::Xlinspace(int resolution)
{
    std::vector<double> steps;
    double left_limit_x = assembler->getBspline_x().getKnotvector()(0);
    double right_limit_x = assembler->getBspline_x().getKnotvector()(assembler->getBspline_x().getKnotvector().getSize() - 1);
    for (int i = (int)(left_limit_x); i <= resolution; i++)
    {
        double i_step = left_limit_x + (double)(i) * ((right_limit_x - left_limit_x) / ((double)(resolution)));
        steps.push_back(i_step);
    }
    return steps;
}

template <class T>
void Poisson<T>::plotSolution(int resolution, std::string filename)
{
    // Create B-spline curve
    std::vector<double> steps = Xlinspace(resolution);

    std::ofstream my_file(filename);
    my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
    my_file << "zone t= " << "\"1\"" << ",i=" << resolution + 1 << ",j=" << resolution + 1 << "\n";

    for (int i = 0; i < steps.size(); i++)
    {
        int span = assembler->XspanOfValueInKnotVector(steps[i]);
        std::vector<double> bVal = assembler->getBspline_x().evaluateAtPoint(steps[i]).first;

        double coord_x = 0.0, coord_y = 0.0, coord_z = 0.0;
        for (int kk = 0; kk < bVal.size(); kk++)
        {
            coord_x += bVal[kk] * assembler->getControlPoints()[span - assembler->getBspline_x().getDegree() + kk][0];
            coord_y += bVal[kk] * assembler->getControlPoints()[span - assembler->getBspline_x().getDegree() + kk][1];
            coord_z += bVal[kk] * solution[span - assembler->getBspline_x().getDegree() + kk]; // interpolation
        }
        my_file << coord_x << " " << coord_z << "\n";
    }
    my_file.close();
}

template <class T>
Solver *Poisson<T>::getSolver()
{
    return solver;
}

template <class T>
std::vector<double> &Poisson<T>::getSolution()
{
    return solution;
}

template <class T>
T *Poisson<T>::getAssembler()
{
    return assembler;
}

template <class T>
void Poisson<T>::setSolution(std::vector<double> &&newSolution)
{ 
    solution = newSolution;
    expandSolutionOnBoundary();
}

template <class T>
void Poisson<T>::setSolution(std::vector<double> &newSolution)
{
    solution = newSolution;
    expandSolutionOnBoundary();
}

template <class T>
void Poisson<T>::updateRhs(std::vector<double> &b)
{
    solver->setLeftAndRightHandSides(assembler->getSystemMatrix(), b);
}