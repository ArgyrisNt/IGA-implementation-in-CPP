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
                    newSolution.push_back(assembler->boundaryConditions->west.second);
                }
                else if (assembler->boundaryBasisFunctions[position].second == 2)
                {
                    newSolution.push_back(assembler->boundaryConditions->east.second);
                }
                else if (assembler->boundaryBasisFunctions[position].second == 3)
                {
                    newSolution.push_back(assembler->boundaryConditions->south.second);
                }
                else if (assembler->boundaryBasisFunctions[position].second == 4)
                {
                    newSolution.push_back(assembler->boundaryConditions->north.second);
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
Solver* Poisson<T>::getSolver()
{
    return solver;
}

template <class T>
std::vector<double> &Poisson<T>::getSolution()
{
    return solution;
}

template <class T>
T* Poisson<T>::getAssembler()
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