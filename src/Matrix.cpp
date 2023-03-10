#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include "..\include\Matrix.h"


template <class T>
Matrix<T>::Matrix()
{
	numberOfRows = 0;
	numberOfColumns = 0;
	values = new T *[numberOfRows];
	for (int i = 0; i < numberOfRows; i++)
	{
		values[i] = new T[numberOfColumns];
	}
}

template <class T>
Matrix<T>::Matrix(int rows, int columns, T value)
{
	numberOfRows = rows;
	numberOfColumns = columns;
	values = new T *[numberOfRows];
	for (int i = 0; i < numberOfRows; i++)
	{
		values[i] = new T[numberOfColumns];
		for (int j = 0; j < numberOfColumns; j++)
		{
			values[i][j] = value;
		}
	}
}

template <class T>
Matrix<T>::Matrix(const Matrix &matrix)
{
	numberOfRows = matrix.numberOfRows;
	numberOfColumns = matrix.numberOfColumns;
	values = matrix.values;
}

template <class T>
Matrix<T>::~Matrix() {}

template <class T>
Matrix<T> &Matrix<T>::operator=(Matrix matrix)
{
	numberOfRows = matrix.numberOfRows;
	numberOfColumns = matrix.numberOfColumns;
	values = matrix.values;

	return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator+(Matrix matrix)
{
	assert(matrix.getNumberOfRows() == getNumberOfRows());
	assert(matrix.getNumberOfColumns() == getNumberOfColumns());
	Matrix<T> result(numberOfRows, numberOfColumns);
	for (int i = 0; i < getNumberOfRows(); i++)
	{
		for (int j = 0; j < getNumberOfColumns(); j++)
		{
			result.setValue(i, j, values[i][j] + matrix.values[i][j]);
		}
	}

	return result;
}

template <class T>
Matrix<T> Matrix<T>::operator-(Matrix &matrix)
{
	assert(matrix.getNumberOfRows() == getNumberOfRows());
	assert(matrix.getNumberOfColumns() == getNumberOfColumns());
	Matrix<T> result(numberOfRows, numberOfColumns);
	for (int i = 0; i < getNumberOfRows(); i++)
	{
		for (int j = 0; j < getNumberOfColumns(); j++)
		{
			result.setValue(i, j, values[i][j] - matrix.values[i][j]);
		}
	}
	return result;
}

template <class T>
Matrix<T> Matrix<T>::operator*(Matrix &matrix)
{
	assert(matrix.getNumberOfRows() == getNumberOfColumns());
	Matrix result(getNumberOfRows(), matrix.getNumberOfColumns());
	for (int i = 0; i < getNumberOfRows(); ++i)
	{
		for (int j = 0; j < matrix.getNumberOfColumns(); ++j)
		{
			for (int k = 0; k < getNumberOfColumns(); ++k)
			{
				result.values[i][j] += values[i][k] * matrix.values[k][j];
			}
		}
	}
	return result;
}

template <class T>
Matrix<T> Matrix<T>::operator*(const T &constant)
{
	Matrix<T> result(numberOfRows, numberOfColumns);
	for (int i = 0; i < getNumberOfRows(); i++)
	{
		for (int j = 0; j < getNumberOfColumns(); j++)
		{
			result.setValue(i, j, values[i][j] * constant);
		}
	}

	return result;
}

template <class T>
std::vector<T> Matrix<T>::operator*(std::vector<T> &matrixVector)
{
	assert(matrixVector.size() == getNumberOfColumns());
	std::vector<T> result(getNumberOfRows(), 0.0);
	for (int i = 0; i < getNumberOfRows(); i++)
	{
		for (int j = 0; j < getNumberOfColumns(); j++)
		{
			result[i] += values[i][j] * matrixVector[j];
		}
	}

	return result;
}

template <class T>
T Matrix<T>::operator()(int row, int column)
{
	assert(row < getNumberOfRows());
	assert(column < getNumberOfColumns());
	return values[row][column];
}

template <class T>
void Matrix<T>::setValue(int row, int column, T value)
{
	assert(row < getNumberOfRows());
	assert(column < getNumberOfColumns());
	values[row][column] = value;
}

template <typename T>
double norm(std::vector<T> &v)
{
	double result = 0.0;
	for (size_t i = 0; i < v.size(); i++)
	{
		result += v[i] * v[i];
	}
	result = sqrt(result);

	return result;
}

template <class T>
T Matrix<T>::determinant()
{
	switch (getNumberOfRows())
	{
	case 1:
		return values[0][0];
		break;
	case 2:
		return (values[0][0] * values[1][1] - values[0][1] * values[1][0]);
		break;
	default:
		throw std::invalid_argument("Invalid dimensions");
		break;
	}
}

template <class T>
void Matrix<T>::print()
{
	for (int i = 0; i < getNumberOfRows(); i++)
	{
		for (int j = 0; j < getNumberOfColumns(); j++)
		{
			std::cout << values[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

template <class T>
Matrix<T> Matrix<T>::transpose()
{
	Matrix result(getNumberOfRows(), getNumberOfColumns());
	for (int i = 0; i < getNumberOfRows(); i++)
	{
		for (int j = 0; j < getNumberOfColumns(); j++)
		{
			result.values[i][j] = values[j][i];
		}
	}

	return result;
}

template <class T>
Matrix<T> Matrix<T>::inverse()
{
	if (determinant() > -1e-7 && determinant() < 1e-7) throw std::invalid_argument("Non invertible matrix");

	Matrix result(getNumberOfRows(), getNumberOfColumns());
	switch (getNumberOfRows())
	{
	case 1:
		result.values[0][0] = 1.0 / values[0][0];
		break;
	case 2:
		result.values[0][0] = (1 / determinant()) * values[1][1];
		result.values[1][1] = (1 / determinant()) * values[0][0];
		result.values[0][1] = -(1 / determinant()) * values[0][1];
		result.values[1][0] = -(1 / determinant()) * values[1][0];
		break;
	default:
		throw std::invalid_argument("Invalid dimensions");
		break;
	}

	return result;
}

template <class T>
std::vector<T> Matrix<T>::Jacobi_iterator(std::vector<T> rightHandSide, int numberOfIterations)
{
	double threshold = 1e-7;
	assert(getNumberOfRows() == rightHandSide.size());
	assert(getNumberOfColumns() == rightHandSide.size());
	std::cout << "\nSolving with Jacobi iterative method. . ." << "\n";
	size_t n = getNumberOfRows();
	std::vector<T> solution(n, 0.0), residual(n), estimation(n);
	int iteration = 0;
	for (iteration = 0; iteration <= numberOfIterations; iteration++)
	{
		std::vector<T> c(n);
		for (size_t i = 0; i < n; i++)
		{
			c[i] = rightHandSide[i];
			for (size_t j = 0; j < n; j++)
			{
				if (j != i)	c[i] -= values[i][j] * solution[j];
			}
		}
		for (size_t i = 0; i < n; i++)
		{
			solution[i] = c[i] / values[i][i];
		}
		estimation = (*this) * solution;
		for (size_t k = 0; k < solution.size(); k++)
		{
			residual[k] = rightHandSide[k] - estimation[k];
		}
		if (norm(residual) < threshold)
		{
			break;
		}
	}

	std::cout << iteration << " iterations\n";
	return solution;
}

template <class T>
std::vector<T> Matrix<T>::GaussSeidel_iterator(std::vector<T> rightHandSide, int numberOfIterations)
{
	double threshold = 1e-7;
	assert(getNumberOfRows() == rightHandSide.size());
	assert(getNumberOfColumns() == rightHandSide.size());
	std::cout << "\nSolving with Gauss Seidel iterative method. . ." << "\n";
	size_t n = getNumberOfRows();
	std::vector<T> solution(n, 0.0), y(n), residual(n), estimation(n);
	int iteration = 0;
	for (iteration = 0; iteration <= numberOfIterations; iteration++)
	{
		for (size_t i = 0; i < n; i++)
		{
			y[i] = rightHandSide[i] / values[i][i];
			for (size_t j = 0; j < n; j++)
			{
				if (j != i)
				{
					y[i] -= (values[i][j] / values[i][i]) * solution[j];
					solution[i] = y[i];
				}
			}
		}
		estimation = (*this) * solution;
		for (size_t k = 0; k < solution.size(); k++)
		{
			residual[k] = rightHandSide[k] - estimation[k];
		}
		if (norm(residual) < threshold)
		{
			break;
		}
	}

	std::cout << iteration << " iterations\n";
	return solution;
}

template <class T>
std::vector<T> Matrix<T>::SOR_iterator(std::vector<T> rightHandSide, int numberOfIterations, double omega)
{
	double threshold = 1e-7;
	assert(getNumberOfRows() == rightHandSide.size());
	assert(getNumberOfColumns() == rightHandSide.size());
	std::cout << "\nSolving with SOR iterative method. . ." << "\n";
	size_t n = getNumberOfRows();
	std::vector<T> solution(n, 0.25), residual(n), estimation(n);
	for (size_t i = 0; i < solution.size(); i++)
	{
		if ((i == 0) || (i == solution.size() - 1)) solution[i] = 0.5;
	}
	int iteration = 0;
	for (iteration = 0; iteration <= numberOfIterations; iteration++)
	{
		for (size_t i = 0; i < n; i++)
		{
			double sigma = 0;
			for (size_t j = 0; j < n; j++)
			{
				if (j != i) sigma += values[i][j] * solution[j];
			}
			solution[i] = (1 - omega) * solution[i] + (omega / values[i][i]) * (rightHandSide[i] - sigma);
		}
		estimation = (*this) * solution;
		for (size_t k = 0; k < solution.size(); k++)
		{
			residual[k] = rightHandSide[k] - estimation[k];
		}
		if (norm(residual) < threshold)
		{
			break;
		}
	}

	std::cout << iteration << " iterations\n";
	return solution;
}

template <class T>
std::vector<T> Matrix<T>::gradient_iterator(std::vector<T> rightHandSide, int numberOfIterations)
{
	std::cout << "\nSolving with Gradient iterative method. . ." << "\n";
	size_t n = getNumberOfRows();
	std::vector<T> solution(n, 0.0), residual(n);
	double nominator = 0.0;
	double denominator = 0.0;
	int iteration = 0;
	for (iteration = 0; iteration < numberOfIterations; iteration++)
	{
		for (size_t i = 0; i < solution.size(); i++)
		{
			residual[i] = rightHandSide[i] - ((*this) * solution)[i];
		}
		for (size_t i = 0; i < solution.size(); i++)
		{
			nominator += residual[i] * residual[i];
			denominator += residual[i] * ((*this) * residual)[i];
		}		
		for (size_t i = 0; i < solution.size(); i++)
		{
			solution[i] = solution[i] + (nominator / denominator) * residual[i];
		}
		if (norm(residual) < 1e-7)
		{
			break;
		}
	}

	std::cout << iteration << " iterations\n";
	return solution;
}

template <typename T>
std::vector<T> Matrix<T>::conjugate_gradient_iterator(std::vector<T> &rightHandSide, int numberOfIterations)
{
	std::cout << "\nSolving with Conjugate Gradient iterative method. . ." << "\n";
	size_t n = getNumberOfRows();
	std::vector<T> solution(n, 0.0), residual(n), t(n);
	double nominator = 0, denominator = 0, nom = 0, denom = 0;
	int iteration = 0;
	for (iteration = 0; iteration < numberOfIterations; iteration++)
	{
		for (size_t i = 0; i < solution.size(); i++)
		{
			residual[i] = rightHandSide[i] - ((*this) * solution)[i];
		}
		if (iteration == 0) t = residual;
		for (size_t i = 0; i < solution.size(); i++)
		{
			nominator += residual[i] * residual[i];
			denominator += t[i] * ((*this) * t)[i];
		}
		for (size_t i = 0; i < solution.size(); i++)
		{
			solution[i] = solution[i] + (nominator / denominator) * t[i];
			residual[i] = residual[i] - (nominator / denominator) * ((*this) * t)[i];
		}
		for (size_t i = 0; i < solution.size(); i++)
		{
			nom += residual[i] * ((*this) * t)[i];
			denom += t[i] * ((*this) * t)[i];
		}
		for (size_t i = 0; i < solution.size(); i++)
		{
			t[i] = residual[i] - (nom / denom) * t[i];
		}
		if (norm(residual) < 1e-7)
		{
			break;
		}
	}

	std::cout << iteration << " iterations\n";
	return solution;
}

template <class T>
std::vector<Matrix<T>> Matrix<T>::LU_factorization()
{
	assert(getNumberOfRows() == getNumberOfColumns());
	std::cout << "\nSolving with LU decomposition method. . ." << "\n";
	size_t n = getNumberOfRows();
	Matrix L(n, n);
	Matrix U(n, n);

	// Doolittle algorithm
	for (size_t i = 0; i < n; i++)
	{
		// matrix U
		for (size_t k = i; k < n; k++)
		{
			T sum = 0.0;
			for (size_t j = 0; j < i; j++)
			{
				sum += L(i, j) * U(j, k);
			}
			U.setValue(i, k, values[i][k] - sum);
		}

		// matrix L
		for (size_t k = i; k < n; k++)
		{
			if (i == k) L.setValue(i, i, 1);
			else
			{
				T sum = 0.0;
				for (size_t j = 0; j < i; j++)
				{
					sum += L(k, j) * U(j, i);
				}
				L.setValue(k, i, (values[k][i] - sum) / U(i, i));
			}
		}
	}
	std::vector<Matrix> result(2, Matrix(n, n));
	result[0] = L;
	result[1] = U;

	return result;
}

template <class T>
std::vector<T> Matrix<T>::forward_Euler(std::vector<T> rightHandSide)
{
	size_t n = getNumberOfRows();
	std::vector<T> solution(n);
	for (size_t j = 0; j < n; j++)
	{
		solution[0] = rightHandSide[0] / values[0][0];
		for (size_t i = 1; i < n; i++)
		{
			T sum = 0;
			for (size_t k = 0; k <= i - 1; k++)
			{
				sum += values[i][k] * solution[k];
			}
			solution[i] = (rightHandSide[i] - sum) / values[i][i];
		}
	}
	return solution;
}

template <class T>
std::vector<T> Matrix<T>::backward_Euler(std::vector<T> rightHandSide)
{
	int n = getNumberOfRows();
	std::vector<T> solution(n);

	for (int j = 0; j < n; j++)
	{
		solution[n - 1] = rightHandSide[n - 1] / values[n - 1][n - 1]; // OK
		for (int i = n - 2; i >= 0; i--)
		{
			T sum = 0.0;
			for (int k = i + 1; k < n; k++)
			{
				sum += values[i][k] * solution[k];
			}
			solution[i] = (rightHandSide[i] - sum) / values[i][i];
		}
	}
	return solution;
}

template <class T>
std::vector<Matrix<T>> Matrix<T>::QR_factorization()
{
	assert(getNumberOfRows() == getNumberOfColumns());
	std::cout << "\nSolving with QR decomposition method. . ." << "\n";
	size_t n = getNumberOfRows();
	Matrix Q(n, n);
	Matrix R(n, n);

	int counter = 0;
	Matrix ATranspose = transpose();
	Matrix QTranspose = Q.transpose();
	for (int j = 0; j < n; j++)
	{
		std::vector<T> e;
		std::vector<T> a;
		for (int column = 0; column < n; column++)
		{
			a.push_back(ATranspose.values[j][column]);
		}
		std::vector<T> u(a);

		for (int i = 0; i < counter; i++)
		{
			T temp = 0.0;
			for (int k = 0; k < n; k++)
			{
				temp += Q(k,i) * a[k];
			}

			std::vector<T> proj;
			for (int k = 0; k < n; k++)
			{
				proj.push_back(Q(k,i) * temp);
			}

			for (int k = 0; k < n; k++)
			{
				u[k] = u[k] - proj[k];
			}
		}

		T result = 0;
		for (size_t i = 0; i < u.size(); i++)
		{
			result += u[i] * u[i];
		}
		result = sqrt(result);

		for (int k = 0; k < n; k++)
		{
			e.push_back(u[k]/result);
		}

		for (int k = 0; k < n; k++)
		{
			Q.setValue(k,counter,e[k]);
		}

		counter++;
	}

	QTranspose = Q.transpose();
	R = QTranspose * (*this);

	std::vector<Matrix> result(2, Matrix(n, n));
	result[0] = Q;
	result[1] = R;

	return result;
}
