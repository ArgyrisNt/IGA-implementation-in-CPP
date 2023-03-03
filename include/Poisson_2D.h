#pragma once

#include <iostream>
#include "..\include\Solver.h"
#include "..\include\Assembler_2D.h"
#include "..\include\Poisson.h"

class Poisson_2D : public Poisson<Assembler_2D>
{
public:
    // Constructor
    Poisson_2D(Assembler_2D& ass, Solver& solv) : Poisson<Assembler_2D>(ass, solv) {}

    // Destructor
    ~Poisson_2D() {}

    // Member function
	void plotSol(int resolution);
};
