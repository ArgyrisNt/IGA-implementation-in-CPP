#include "..\include\Solver.h"

std::string name1("LU");
Solver Solver::LU(name1);

std::string name2 = "QR";
Solver Solver::QR(name2);

std::string name3 = "GaussSeidel";
Solver Solver::GaussSeidel(name3);

std::string name4 = "Jacobi";
Solver Solver::Jacobi(name4);

std::string name5 = "SOR";
Solver Solver::SOR(name5);

std::string name6 = "Gradient";
Solver Solver::Gradient(name6);

std::string name7 = "ConjugateGradient";
Solver Solver::ConjugateGradient(name7);

void Solver::LUsolve(Matrix<double> &A, std::vector<double> &b)
{
    std::vector<Matrix<double>> LU = A.LU_factorization();
    Matrix<double> L = LU[0];
    Matrix<double> U = LU[1];
    std::vector<double> y;
    y = L.forward_Euler(b);
    solution = U.backward_Euler(y);
}

void Solver::QRsolve(Matrix<double> &A, std::vector<double> &b)
{
    std::vector<Matrix<double>> QR = A.QR_factorization();
    std::vector<double> rhs = QR[0].transpose() * b;
    solution = QR[1].backward_Euler(rhs);
}