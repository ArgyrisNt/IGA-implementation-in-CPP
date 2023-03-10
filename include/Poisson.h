#pragma once

#include <iostream>
#include "..\include\Solver.h"

template<class T>
class Poisson
{
public:
    // Constructor
    Poisson(T &newAssembler, Solver &newSolver) : assembler(&newAssembler), solver(&newSolver){};

    // Destructor
    virtual ~Poisson() {}

    // Member functions
    virtual void expandSolutionOnBoundary();
    virtual void plotSolution(int resolution);

    // Member setter functions
    void setSolution(std::vector<double> &newSolution) { solution = newSolution; }

    // Member getter functions
    std::vector<double>& getSolution() { return solution; }
    T* getAssembler() { return assembler; }
    Solver* getSolver() { return solver; }

protected:
    // Member variables
    std::vector<double> solution;
    T *assembler;
    Solver *solver;
};

inline std::ostream& operator << (std::ostream& os, std::vector<double>& vec)
{
    os << "Solution is: " << std::endl;
    for (auto value : vec)
    {
        os << value << std::endl;
    }

    return os;
}