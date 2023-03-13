#ifndef H_POISSON
#define H_POISSON

#include <iostream>
#include "..\include\Solver.h"

template<class T>
class Poisson
{
public:
    Poisson(T &newAssembler, Solver &newSolver);

    virtual ~Poisson() {}

    virtual void expandSolutionOnBoundary();
    virtual void plotSolution(int resolution, std::string filename);

    void setSolution(std::vector<double> &&newSolution);
    void setSolution(std::vector<double> &newSolution);

    std::vector<double>& getSolution();
    T* getAssembler();
    Solver* getSolver();

    void updateRhs(std::vector<double>& b);

    std::vector<double> Xlinspace(int resolution);
    
protected:
    std::vector<double> solution;
    T *assembler;
    Solver *solver;
};

inline std::ostream& operator << (std::ostream& os, std::vector<double>& vec)
{
    for (auto value : vec)
    {
        os << value << std::endl;
    }

    return os;
}

#include "..\src\Poisson.cpp"

#endif