#pragma once

#include <iostream>
#include "..\include\Solver.h"

template<class T>
class Poisson
{
public:
    // Constructor
    Poisson(T& ass, Solver& solv) : assembler(&ass), solver(&solv) {};
    
    // Destructor
    virtual ~Poisson() {}

    // Member functions
    virtual void constructSol();
    virtual void plotSol(std::string filename1, std::string filename2);

    // Member setter functions
    void setSolution(std::vector<double>& _solution) { solution = _solution; }

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

inline std::ostream& operator << (std::ostream& os, std::vector<double>& obj)
{
    os << "Solution is: " << std::endl;
    for (auto el: obj)
    {
        os << el << std::endl;
    }

    return os;
}