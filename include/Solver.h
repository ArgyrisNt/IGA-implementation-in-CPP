#pragma once

#include <iostream>
#include <vector>
#include "..\include\Matrix.h"

class Solver
{
public:
    // Common solvers
    static Solver LU;
    static Solver QR;
    static Solver GaussSeidel;
    static Solver Jacobi;
    static Solver SOR;
    static Solver Gradient;
    static Solver ConjugateGradient;

    // Constructor
    Solver(const std::string& name) : mode(name) {}

    // Destructor
    ~Solver() {}

    // Member functions
    std::vector<double>& solve(Matrix<double>& A, std::vector<double>& b, int numberOfIterations = 50, double omega = 1.03);

    // Member getter functions
    const std::string getMode() { return mode; }

private:
    void LUsolve(Matrix<double> &leftHandSide, std::vector<double> &rightHandSide);
    void QRsolve(Matrix<double> &leftHandSide, std::vector<double> &rightHandSide);


    // Member variables
    std::vector<double> solution;
    const std::string mode;
};