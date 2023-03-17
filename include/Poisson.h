#ifndef H_POISSON
#define H_POISSON

#include <iostream>
#include "..\include\Solver.h"

template<class T>
class Poisson
{
public:
    Poisson(T &newAssembler, Solver &newSolver);

    ~Poisson() {}

    void expandSolutionOnBoundary();

    void setSolution(std::vector<double> &&newSolution);
    void setSolution(std::vector<double> &newSolution);

    std::vector<double>& getSolution();
    T* getAssembler();
    Solver* getSolver();

    void updateRhs(std::vector<double>& b);
    
protected:
    std::vector<double> solution;
    T* assembler;
    Solver* solver;
};

#include "..\src\Poisson.cpp"

#endif