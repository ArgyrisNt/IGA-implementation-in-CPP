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
	for (int i = 0; i < numberOfRows; ++i)
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
	for (int i = 0; i < numberOfRows; ++i)
	{
		values[i] = new T[numberOfColumns];
		for (int j = 0; j < numberOfColumns; ++j)
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
Matrix<T> &Matrix<T>::operator=(const Matrix &matrix)
{
	numberOfRows = matrix.numberOfRows;
	numberOfColumns = matrix.numberOfColumns;
	values = matrix.values;

	return *this;
}


template <class T>
void Matrix<T>::setValue(int row, int column, T value)
{
	assert(row < getNumberOfRows());
	assert(column < getNumberOfColumns());
	values[row][column] = value;
}


template <class T>
T Matrix<T>::determinant() const
{
	assert(getNumberOfRows() == getNumberOfColumns());
	int n = getNumberOfRows();
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
void Matrix<T>::print() const
{
	for (int i = 0; i < getNumberOfRows(); ++i)
	{
		for (int j = 0; j < getNumberOfColumns(); ++j)
		{
			std::cout << values[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

template <class T>
Matrix<T> Matrix<T>::transpose() const
{
	Matrix result(getNumberOfRows(), getNumberOfColumns());
	for (int i = 0; i < getNumberOfRows(); ++i)
	{
		for (int j = 0; j < getNumberOfColumns(); ++j)
		{
			result.values[i][j] = values[j][i];
		}
	}

	return result;
}

template <class T>
Matrix<T> Matrix<T>::inverse() const
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
std::vector<Matrix<T>> Matrix<T>::LU_factorization()
{
	assert(getNumberOfRows() == getNumberOfColumns());
	std::cout << "\nSolving with LU decomposition method. . ." << "\n";
	size_t n = getNumberOfRows();
	Matrix L(n, n);
	Matrix U(n, n);

	// Doolittle algorithm
	for (size_t i = 0; i < n; ++i)
	{
		// matrix U
		for (size_t k = i; k < n; ++k)
		{
			T sum = 0.0;
			for (size_t j = 0; j < i; ++j)
			{
				sum += L(i, j) * U(j, k);
			}
			U.setValue(i, k, values[i][k] - sum);
		}

		// matrix L
		for (size_t k = i; k < n; ++k)
		{
			if (i == k) L.setValue(i, i, 1);
			else
			{
				T sum = 0.0;
				for (size_t j = 0; j < i; ++j)
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
std::vector<T> Matrix<T>::forward_Euler(const std::vector<T> &rightHandSide)
{
	size_t n = getNumberOfRows();
	std::vector<T> solution(n);
	for (size_t j = 0; j < n; ++j)
	{
		solution[0] = rightHandSide[0] / values[0][0];
		for (size_t i = 1; i < n; ++i)
		{
			T sum = 0;
			for (size_t k = 0; k <= i - 1; ++k)
			{
				sum += values[i][k] * solution[k];
			}
			solution[i] = (rightHandSide[i] - sum) / values[i][i];
		}
	}
	return solution;
}

template <class T>
std::vector<T> Matrix<T>::backward_Euler(const std::vector<T> &rightHandSide)
{
	int n = getNumberOfRows();
	std::vector<T> solution(n);

	for (int j = 0; j < n; ++j)
	{
		solution[n - 1] = rightHandSide[n - 1] / values[n - 1][n - 1]; // OK
		for (int i = n - 2; i >= 0; i--)
		{
			T sum = 0.0;
			for (int k = i + 1; k < n; ++k)
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
	for (int j = 0; j < n; ++j)
	{
		std::vector<T> e;
		std::vector<T> a;
		for (int column = 0; column < n; ++column)
		{
			a.push_back(ATranspose.values[j][column]);
		}
		std::vector<T> u(a);

		for (int i = 0; i < counter; ++i)
		{
			T temp = 0.0;
			for (int k = 0; k < n; ++k)
			{
				temp += Q(k,i) * a[k];
			}

			std::vector<T> proj;
			for (int k = 0; k < n; ++k)
			{
				proj.push_back(Q(k,i) * temp);
			}

			for (int k = 0; k < n; ++k)
			{
				u[k] = u[k] - proj[k];
			}
		}

		T result = 0;
		for (size_t i = 0; i < u.size(); ++i)
		{
			result += u[i] * u[i];
		}
		result = sqrt(result);

		for (int k = 0; k < n; ++k)
		{
			e.push_back(u[k]/result);
		}

		for (int k = 0; k < n; ++k)
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
