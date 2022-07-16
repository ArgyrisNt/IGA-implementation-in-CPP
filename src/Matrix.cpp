#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include "C:\Users\argir\OneDrive\Desktop\IGA_new\include\Matrix.h"

Matrix::Matrix()
{
	rows = 0;
	cols = 0;
	detVal = 0;
	mat = new double* [rows];
	for (int i = 0; i < rows; i++)
	{
		mat[i] = new double[cols];
	}
}

Matrix::Matrix(int my_rows, int my_cols)
{
	rows = my_rows;
	cols = my_cols;
	detVal = 0;
	mat = new double* [rows];
	for (int i = 0; i < rows; i++)
	{
		mat[i] = new double[cols];
	}
}

Matrix::Matrix(int my_rows, int my_cols, double val)
{
	rows = my_rows;
	cols = my_cols;
	detVal = 0;
	mat = new double* [rows];
	for (int i = 0; i < rows; i++)
	{
		mat[i] = new double[cols];
	}
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			mat[i][j] = val;
		}
	}
}

Matrix::Matrix(const Matrix& m)
{
	rows = m.rows;
	cols = m.cols;
	detVal = m.detVal;
	mat = m.mat;
}

Matrix::~Matrix() {}

Matrix& Matrix::operator=(Matrix m)
{
	rows = m.rows;
	cols = m.cols;
	detVal = m.detVal;
	mat = m.mat;

	return *this;
}

Matrix Matrix::operator+(Matrix& m)
{
	assert(m.rows == rows);
	assert(m.cols == cols);
	Matrix res(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			res.setValue(i, j, mat[i][j] + m.mat[i][j]);
		}
	}

	return res;
}

Matrix Matrix::operator-(Matrix& m)
{
	assert(m.rows == rows);
	assert(m.cols == cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			mat[i][j] -= m.mat[i][j];
		}
	}
	return *this;
}

Matrix Matrix::operator*(Matrix& m)
{
	assert(m.rows == cols);
	Matrix n(rows, m.cols);
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < m.cols; ++j)
		{
			n.setValue(i, j, 0.0);
		}
	}
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < m.cols; ++j)
		{
			for (int k = 0; k < cols; ++k)
			{
				n.mat[i][j] += mat[i][k] * m.mat[k][j];
			}
		}
	}
	return n;
}

Matrix Matrix::operator*(const double& c)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			mat[i][j] = mat[i][j] * c;
		}
	}

	return *this;
}

std::vector<double> Matrix::operator*(std::vector<double>& vec)
{
	assert(vec.size() == cols);
	std::vector<double> res(rows, 0.0);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			res[i] += mat[i][j] * vec[j];
		}
	}

	return res;
}

double Matrix::operator()(int row, int col)
{
	assert(row < rows);
	assert(col < cols);
	return mat[row][col];
}

void Matrix::setValue(int row, int col, double val)
{
	assert(row < rows);
	assert(col < cols);
	mat[row][col] = val;
}

void Matrix::calcDet()
{
	switch (rows)
	{
	case 1:
		detVal = mat[0][0];
		break;
	case 2:
		detVal = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		break;
	default:
		throw std::invalid_argument("Invalid dimensions");
		break;
	}
}

void Matrix::print()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::cout << mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

Matrix Matrix::transpose()
{
	Matrix result(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			result.mat[i][j] = mat[j][i];
		}
	}
	result.detVal = detVal;

	return result;
}

Matrix Matrix::inverse()
{
	Matrix result(rows, cols);
	switch (rows)
	{
	case 1:
		result.mat[0][0] = 1.0 / mat[0][0];
		break;
	case 2:
		calcDet();
		if (!detVal)
		{
			throw std::invalid_argument("Non invertible matrix");
			break;
		}
		result.detVal = 1 / detVal;
		result.mat[0][0] = result.detVal * mat[1][1];
		result.mat[1][1] = result.detVal * mat[0][0];
		result.mat[0][1] = -result.detVal * mat[0][1];
		result.mat[1][0] = -result.detVal * mat[1][0];
		break;
	default:
		throw std::invalid_argument("Invalid dimensions");
		break;
	}

	return result;
}

double Matrix::norm(std::vector<double>& v)
{
	double result = 0;
	for (size_t i = 0; i < v.size(); i++)
	{
		result += v[i] * v[i];
	}
	result = sqrt(result);

	return result;
}

std::vector<double> Matrix::Jacobi_iterator(std::vector<double> b, int iters)
{
	assert(rows == b.size());
	assert(cols == b.size());
	std::cout << "\nSolving with Jacobi iterative method. . ." << "\n";
	size_t n = rows;
	std::vector<double> x(n), r(n), sec(n);
	for (size_t i = 0; i < x.size(); i++)
	{
		x[i] = 0.0;
	}
	double sum1 = 0;
	double sum2 = 0;
	for (int iter = 0; iter <= iters; iter++)
	{
		std::vector<double> c(n);
		for (size_t i = 0; i < n; i++)
		{
			c[i] = b[i];
			for (size_t j = 0; j < n; j++)
			{
				if (j != i)	c[i] -= mat[i][j] * x[j];
			}
		}
		for (size_t i = 0; i < n; i++)
		{
			x[i] = c[i] / mat[i][i];
		}
	}
	return x;
}

std::vector<Matrix> Matrix::LU_factor() 
{
	assert(rows == cols);
	std::cout << "\nSolving with LU decomposition method. . ." << "\n";
	size_t n = rows;
	Matrix L(n, n);
	Matrix U(n, n);

	// Doolittle algorithm
	for (size_t i = 0; i < n; i++)
	{
		// matrix U
		for (size_t k = i; k < n; k++)
		{
			double sum = 0.0;
			for (size_t j = 0; j < i; j++)
			{
				sum += L(i, j) * U(j, k);
			}
			U.setValue(i, k, mat[i][k] - sum);
		}

		// matrix L
		for (size_t k = i; k < n; k++)
		{
			if (i == k) L.setValue(i, i, 1);
			else
			{
				double sum = 0.0;
				for (size_t j = 0; j < i; j++)
				{
					sum += L(k, j) * U(j, i);
				}
				L.setValue(k, i, (mat[k][i] - sum) / U(i, i));
			}
		}
	}
	std::vector<Matrix> result(2, Matrix(n, n));
	result[0] = L;
	result[1] = U;

	return result;
}

std::vector<double> Matrix::forward_Euler(std::vector<double> b) 
{
	size_t n = rows;
	std::vector<double> y(n);
	for (size_t j = 0; j < n; j++)
	{
		y[0] = b[0] / mat[0][0];
		for (size_t i = 1; i < n; i++)
		{
			double sum = 0;
			for (size_t k = 0; k <= i - 1; k++)
			{
				sum += mat[i][k] * y[k];
			}
			y[i] = (b[i] - sum) / mat[i][i];
		}
	}
	return y;
}

std::vector<double> Matrix::backward_Euler(std::vector<double> b)
{
	int n = rows;
	std::vector<double> x(n);

	for (int j = 0; j < n; j++)
	{
		x[n - 1] = b[n - 1] / mat[n - 1][n - 1]; // OK
		for (int i = n - 2; i >= 0; i--)
		{
			double sum = 0.0;
			for (int k = i + 1; k < n; k++)
			{
				sum += mat[i][k] * x[k];
			}
			x[i] = (b[i] - sum) / mat[i][i];
		}
	}
	return x;
}