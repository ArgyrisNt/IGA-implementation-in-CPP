#include <iostream>
#include "..\include\Poisson.h"

template <class T>
Poisson<T>::Poisson(T &newAssembler, Solver &newSolver)
{
    assembler = std::make_shared<T>(newAssembler);
    solver = std::make_shared<Solver>(newSolver);
    solver->setLeftAndRightHandSides(assembler->getSystemMatrix(), assembler->getRightHandSide());
}



template <class T>
void Poisson<T>::applyInitialCondition(std::vector<double>& initialSolution)
{
    solution = initialSolution;
}

template <class T>
void Poisson<T>::expandSolutionOnBoundary()
{
    std::vector<double> newSolution;
    if (assembler->getBoundarymode() == "Multipliers")
    {
        solution.erase(solution.begin(), solution.begin() + assembler->getBoundaryBasisFunctions().size());
    }
    else
    {
	    int j = 0;
	    for (int i = 0; i < assembler->getNumberOfBasisFunctions(); ++i)
	    {
            auto it = std::find_if(assembler->getBoundaryBasisFunctions().begin(), assembler->getBoundaryBasisFunctions().end(), CompareFirst(i));
            if (it != assembler->getBoundaryBasisFunctions().end())
            {
                int position = it - assembler->getBoundaryBasisFunctions().begin();
                if (assembler->getBoundaryBasisFunctions()[position].second == 1)
                {
                    newSolution.push_back(assembler->getBoundaryConditions().west.second);
                }
                else if (assembler->getBoundaryBasisFunctions()[position].second == 2)
                {
                    newSolution.push_back(assembler->getBoundaryConditions().east.second);
                }
                else if (assembler->getBoundaryBasisFunctions()[position].second == 3)
                {
                    newSolution.push_back(assembler->getBoundaryConditions().south.second);
                }
                else if (assembler->getBoundaryBasisFunctions()[position].second == 4)
                {
                    newSolution.push_back(assembler->getBoundaryConditions().north.second);
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
void Poisson<T>::updateRhs(const std::vector<double> &b)
{
    solver->setLeftAndRightHandSides(assembler->getSystemMatrix(), b);
}

template <class T>
std::shared_ptr<Solver> Poisson<T>::getSolver() const
{
    return solver;
}

template <class T>
std::vector<double> &Poisson<T>::getSolution()
{
    return solution;
}

template <class T>
std::shared_ptr<T> Poisson<T>::getAssembler() const
{
    return assembler;
}



template <class T>
void Poisson<T>::solve(int numberOfIterations, double omega)
{
    solution = solver->solve(numberOfIterations, omega);
    expandSolutionOnBoundary();
}
