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
    std::vector<double>& solve(Matrix<double>& A, std::vector<double>& b, int numberOfIterations = 50, double omega = 1.03)
    {
        if (mode == "LU")
        {
            LUsolve(A, b);
        }
        else if (mode == "QR")
        {
            QRsolve(A, b);
        }
        else if (mode == "GaussSeidel")
        {
            solution = A.GaussSeidel_iterator(b, numberOfIterations);
        }
        else if (mode == "Jacobi")
        {
            solution = A.Jacobi_iterator(b, numberOfIterations);
        }
        else if (mode == "SOR")
        {
            solution = A.SOR_iterator(b, numberOfIterations, omega);
        }
        else if (mode == "Gradient")
        {
            solution = A.gradient_iterator(b, numberOfIterations);
        }
        else if (mode == "ConjugateGradient")
        {
            solution = A.conjugate_gradient_iterator(b, numberOfIterations);      
        }
        else
        {
            std::cout << "Invalid solving method" << std::endl;
		    throw std::invalid_argument("Invalid method");
        }

        return solution;
    }

    // Member getter functions
    const std::string getMode() { return mode; }

private:
    void LUsolve(Matrix<double> &A, std::vector<double> &b);
    void QRsolve(Matrix<double> &A, std::vector<double> &b);

    // Member variables
    std::vector<double> solution;
    const std::string mode;
};