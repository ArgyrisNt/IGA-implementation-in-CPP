#ifndef H_POISSON
#define H_POISSON

#include <iostream>
#include <memory>
#include "..\include\Solver.h"

template<class T>
class Poisson
{
public:
    Poisson(T &newAssembler, Solver &newSolver);

    ~Poisson() {}

    void applyInitialCondition(std::vector<double>&);
    void expandSolutionOnBoundary();
    void solve(int numberOfIterations = 50, double omega = 1.03);

    std::vector<double>& getSolution();
    std::shared_ptr<T> getAssembler() const;
    std::shared_ptr<Solver> getSolver() const;

    void updateRhs(const std::vector<double>& b);
    
private:
    std::vector<double> solution;
    std::shared_ptr<T> assembler;
    std::shared_ptr<Solver> solver;
};

#include "..\src\Poisson.cpp"

#endif