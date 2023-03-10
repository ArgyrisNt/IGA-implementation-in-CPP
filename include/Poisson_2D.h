#pragma once

#include <iostream>
#include "..\include\Solver.h"
#include "..\include\Assembler_2D.h"
#include "..\include\Poisson.h"

class Poisson_2D : public Poisson<Assembler_2D>
{
public:
    // Constructor
    Poisson_2D(Assembler_2D& assembler, Solver& solver) : Poisson<Assembler_2D>(assembler, solver) {}

    // Destructor
    ~Poisson_2D() {}

    // Member function
	void plotSolution(int resolution);
    void writeParameterSpaceToFile(std::string filename);
    void writeTrimmedTrianglesToFile(std::string filename);
};
