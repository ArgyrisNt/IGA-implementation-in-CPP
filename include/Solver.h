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
    std::vector<double>& solve(Matrix<double>& A, std::vector<double>& b, int iters = 10, double omega = 1.03)
    {
        if (mode == "LU")
        {
            std::vector<Matrix<double>> LU = A.LU_factor();
	        Matrix<double> L = LU[0];
            Matrix<double> U = LU[1];
            std::vector<double> y;
	        y = L.forward_Euler(b);
	        solution = U.backward_Euler(y);
            return solution;
        }
        else if (mode == "QR")
        {
            std::vector<Matrix<double>> QR = A.QR_factor();
	        std::vector<double> rhs = QR[0].transpose() * b;
	        std::vector<double> sol = QR[1].backward_Euler(rhs);
            solution = sol;
            return solution;
        }
        else if (mode == "GaussSeidel")
        {
            solution = A.GaussSeidel_iterator(b, iters);
            return solution;
        }
        else if (mode == "Jacobi")
        {
            solution = A.Jacobi_iterator(b, iters);
            return solution;
        }
        else if (mode == "SOR")
        {
            solution = A.SOR_iterator(b, iters, omega);
            return solution;
        }
        else if (mode == "Gradient")
        {
            solution = A.gradient_iterator(b, iters);
            return solution;
        }
        else if (mode == "ConjugateGradient")
        {
            solution = A.conjugate_gradient_iter(b, iters);
            return solution;
        }
        else
        {
            std::cout << "Invalid solving method" << std::endl;
		    throw std::invalid_argument("Invalid method");
        }
    }

    // Member getter functions
    const std::string getMode() { return mode; }

private:
    // Member variables
    const Matrix<double> lhs;
    const std::vector<double> rhs;
    std::vector<double> solution;
    const std::string mode;
};