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