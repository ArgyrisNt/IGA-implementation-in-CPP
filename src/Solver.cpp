#include "..\include\Solver.h"
#include <cassert>

Solver Solver::LU("LU");
Solver Solver::QR("QR");
Solver Solver::GaussSeidel("GaussSeidel");
Solver Solver::Jacobi("Jacobi");
Solver Solver::SOR("SOR");
Solver Solver::Gradient("Gradient");
Solver Solver::ConjugateGradient("ConjugateGradient");

void Solver::setLeftAndRightHandSides(const Matrix<double> &left, const std::vector<double> &right)
{
    leftHandSide = left;
    rightHandSide = right;
}

std::vector<double> Solver::solve(int numberOfIterations, double omega)
{
    assert(leftHandSide.getNumberOfRows() == rightHandSide.size());

    std::vector<double> solution;
    if (mode == "LU") solution = LUsolve();
    else if (mode == "QR") solution = QRsolve();
    else if (mode == "GaussSeidel") solution = GaussSeidel_iterator(numberOfIterations);
    else if (mode == "Jacobi") solution = Jacobi_iterator(numberOfIterations);
    else if (mode == "SOR") solution = SOR_iterator(numberOfIterations, omega);
    else if (mode == "Gradient") solution = gradient_iterator(numberOfIterations);
    else if (mode == "ConjugateGradient") solution = conjugate_gradient_iterator(numberOfIterations);
    else
    {
        std::cout << "Invalid solving method" << std::endl;
        throw std::invalid_argument("Invalid method");
    }

    return solution;
}

std::vector<double> Solver::LUsolve()
{
    std::vector<Matrix<double>> LU = leftHandSide.LU_factorization();
    Matrix<double> L = LU[0];
    Matrix<double> U = LU[1];
    std::vector<double> y = forward_Euler(L, rightHandSide);
    std::vector<double> solution = backward_Euler(U, y);

    return solution;
}

std::vector<double> Solver::QRsolve()
{
    std::vector<Matrix<double>> QR = leftHandSide.QR_factorization();
    std::vector<double> rhs = QR[0].transpose() * rightHandSide;
    std::vector<double> solution = backward_Euler(QR[1], rhs);

    return solution;
}

std::vector<double> Solver::Jacobi_iterator(int numberOfIterations)
{
    double threshold = 1e-7;
    size_t n = leftHandSide.getNumberOfRows();
    std::cout << "\nSolving with Jacobi iterative method. . ." << "\n";
    std::vector<double> solution(n, 0.0), residual(n), estimation(n);
    int iteration = 0;
    for (iteration = 0; iteration <= numberOfIterations; ++iteration)
    {
        std::vector<double> c(n);
        for (size_t i = 0; i < n; ++i)
        {
            c[i] = rightHandSide[i];
            for (size_t j = 0; j < n; ++j)
            {
                if (j == i) continue;
                c[i] -= leftHandSide(i,j) * solution[j];
            }
        }
        for (size_t i = 0; i < n; ++i)
            solution[i] = c[i] / leftHandSide(i,i);

        estimation = leftHandSide * solution;
        for (size_t k = 0; k < solution.size(); ++k)
            residual[k] = rightHandSide[k] - estimation[k];

        if (norm(residual) < threshold) break;
    }

    std::cout << iteration - 1 << " iterations\n";
    return solution;
}

std::vector<double> Solver::GaussSeidel_iterator(int numberOfIterations)
{
    double threshold = 1e-7;
    size_t n = leftHandSide.getNumberOfRows();
    std::cout << "\nSolving with Gauss Seidel iterative method. . ." << "\n";
    std::vector<double> solution(n, 0.0), y(n), residual(n), estimation(n);
    int iteration = 0;
    for (iteration = 0; iteration <= numberOfIterations; ++iteration)
    {
        for (size_t i = 0; i < n; ++i)
        {
            y[i] = rightHandSide[i] / leftHandSide(i,i);
            for (size_t j = 0; j < n; ++j)
            {
                if (j == i) continue;
                y[i] -= (leftHandSide(i,j) / leftHandSide(i,i)) * solution[j];
                solution[i] = y[i];
            }
        }

        estimation = leftHandSide * solution;
        for (size_t k = 0; k < solution.size(); ++k)
            residual[k] = rightHandSide[k] - estimation[k];

        if (norm(residual) < threshold) break;
    }

    std::cout << iteration - 1 << " iterations\n";
    return solution;
}

std::vector<double> Solver::SOR_iterator(int numberOfIterations, double omega)
{
    double threshold = 1e-7;
    size_t n = leftHandSide.getNumberOfRows();
    std::cout << "\nSolving with SOR iterative method. . ." << "\n";
    std::vector<double> solution(n, 0.25), residual(n), estimation(n);
    for (size_t i = 0; i < solution.size(); ++i)
    {
        if ((i == 0) || (i == solution.size() - 1))
            solution[i] = 0.5;
    }
    int iteration = 0;
    for (iteration = 0; iteration <= numberOfIterations; ++iteration)
    {
        for (size_t i = 0; i < n; ++i)
        {
            double sigma = 0;
            for (size_t j = 0; j < n; ++j)
            {
                if (j == i) continue;
                sigma += leftHandSide(i,j) * solution[j];
            }
            solution[i] = (1 - omega) * solution[i] + (omega / leftHandSide(i,i)) * (rightHandSide[i] - sigma);
        }

        estimation = leftHandSide * solution;
        for (size_t k = 0; k < solution.size(); ++k)
            residual[k] = rightHandSide[k] - estimation[k];

        if (norm(residual) < threshold) break;
    }

    std::cout << iteration - 1 << " iterations\n";
    return solution;
}

std::vector<double> Solver::gradient_iterator(int numberOfIterations)
{
    std::cout << "\nSolving with Gradient iterative method. . ." << "\n";
    size_t n = leftHandSide.getNumberOfRows();
    std::vector<double> solution(n, 0.0), residual(n);
    double nominator = 0.0;
    double denominator = 0.0;
    int iteration = 0;
    for (iteration = 0; iteration < numberOfIterations; ++iteration)
    {
        for (size_t i = 0; i < solution.size(); ++i)
            residual[i] = rightHandSide[i] - (leftHandSide * solution)[i];

        for (size_t i = 0; i < solution.size(); ++i)
        {
            nominator += residual[i] * residual[i];
            denominator += residual[i] * (leftHandSide * residual)[i];
        }

        for (size_t i = 0; i < solution.size(); ++i)
            solution[i] = solution[i] + (nominator / denominator) * residual[i];

        if (norm(residual) < 1e-7) break;
    }

    std::cout << iteration - 1 << " iterations\n";
    return solution;
}

std::vector<double> Solver::conjugate_gradient_iterator(int numberOfIterations)
{
    size_t n = leftHandSide.getNumberOfRows();
    std::cout << "\nSolving with Conjugate Gradient iterative method. . ." << "\n";
    std::vector<double> solution(n, 0.0), residual(n), t(n);
    double nominator = 0, denominator = 0, nom = 0, denom = 0;
    int iteration = 0;
    for (iteration = 0; iteration < numberOfIterations; ++iteration)
    {
        for (size_t i = 0; i < solution.size(); ++i)
            residual[i] = rightHandSide[i] - (leftHandSide * solution)[i];

        if (iteration == 0) t = residual;
        for (size_t i = 0; i < solution.size(); ++i)
        {
            nominator += residual[i] * residual[i];
            denominator += t[i] * (leftHandSide * t)[i];
        }
        for (size_t i = 0; i < solution.size(); ++i)
        {
            solution[i] = solution[i] + (nominator / denominator) * t[i];
            residual[i] = residual[i] - (nominator / denominator) * (leftHandSide * t)[i];
        }
        for (size_t i = 0; i < solution.size(); ++i)
        {
            nom += residual[i] * (leftHandSide * t)[i];
            denom += t[i] * (leftHandSide * t)[i];
        }

        for (size_t i = 0; i < solution.size(); ++i)
            t[i] = residual[i] - (nom / denom) * t[i];

        if (norm(residual) < 1e-7) break;
    }

    std::cout << iteration - 1 << " iterations\n";
    return solution;
}

std::vector<double> Solver::forward_Euler(const Matrix<double>& A, const std::vector<double> &b)
{
    size_t n = A.getNumberOfRows();
    assert(n == b.size());

    std::vector<double> solution(n);
    for (size_t j = 0; j < n; ++j)
    {
        solution[0] = b[0] / A(0,0);
        for (size_t i = 1; i < n; ++i)
        {
            double sum = 0;
            for (size_t k = 0; k <= i - 1; ++k)
                sum += A(i,k) * solution[k];
            solution[i] = (b[i] - sum) / A(i,i);
        }
    }
    return solution;
}

std::vector<double> Solver::backward_Euler(const Matrix<double> &A, const std::vector<double> &b)
{
    int n = A.getNumberOfRows();
    assert(n == b.size());
    std::vector<double> solution(n);
    for (int j = 0; j < n; ++j)
    {
        solution[n - 1] = b[n - 1] / A(n - 1,n - 1); // OK
        for (int i = n - 2; i >= 0; i--)
        {
            double sum = 0.0;
            for (int k = i + 1; k < n; ++k)
                sum += A(i,k) * solution[k];
            solution[i] = (b[i] - sum) / A(i,i);
        }
    }
    return solution;
}