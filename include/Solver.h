#ifndef H_SOLVER
#define H_SOLVER

#include <iostream>
#include <vector>
#include "..\include\Matrix.h"

class Solver
{
public:
    static Solver LU, QR, GaussSeidel, Jacobi, SOR, Gradient, ConjugateGradient;

    Solver(const std::string& name) : mode(name) {}

    ~Solver() {}

    std::vector<double> solve(int numberOfIterations = 50, double omega = 1.03);

    void setLeftAndRightHandSides(const Matrix<double>& left, const std::vector<double>& right);
    
private:
    std::vector<double> LUsolve();
    std::vector<double> QRsolve();
    std::vector<double> Jacobi_iterator(int iters);
    std::vector<double> GaussSeidel_iterator(int iters);
    std::vector<double> SOR_iterator(int iters, double omega);
    std::vector<double> gradient_iterator(int iters);
    std::vector<double> conjugate_gradient_iterator(int iters);

    Matrix<double> leftHandSide;
    std::vector<double> rightHandSide;
    const std::string mode;
};

#include "..\src\Solver.cpp"

#endif