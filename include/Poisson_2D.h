#ifndef H_POISSON_2D
#define H_POISSON_2D

#include <iostream>
#include "..\include\Solver.h"
#include "..\include\Assembler_2D.h"
#include "..\include\Poisson.h"

class Poisson_2D : public Poisson<Assembler_2D>
{
public:
    Poisson_2D(Assembler_2D& assembler, Solver& solver) : Poisson<Assembler_2D>(assembler, solver) {}

    ~Poisson_2D() {}

	void plotSolution(int resolution);
    void writeTrimmedTrianglesToFile(std::string filename);

    std::vector<double> Ylinspace(int resolution);
};

#include "..\src\Poisson_2D.cpp"

#endif