#ifndef H_POISSON
#define H_POISSON

#include <iostream>
#include <memory>
#include "Solver.h"

template<class T>
class Poisson
{
public:
    Poisson(T &newAssembler, Solver &newSolver);

    ~Poisson() {}

    void solve(int numberOfIterations = 50, double omega = 1.03);

    void applyInitialCondition(std::vector<double> &initialSolution)
    { solution = initialSolution; }

    void updateRhs(const std::vector<double> &b) 
    { solver->setLeftAndRightHandSides(assembler->getSystemMatrix(), b); }

    const std::vector<double>& getSolution() const 
    { return solution; }
    
private:
    void expandSolutionOnBoundary();

    std::vector<double> solution;
    std::shared_ptr<T> assembler;
    std::shared_ptr<Solver> solver;
};

#include "..\src\Poisson.cpp"

#endif