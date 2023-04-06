#ifndef H_MATRIX
#define H_MATRIX

#include <iostream>
#include "Utilities.h"

using namespace Utils;

template<class T>
class Matrix
{
public:
	Matrix(); 
	Matrix(int my_rows, int my_cols, T val = 0.0);
	Matrix(const Matrix&);

	~Matrix();

	Matrix& operator=(const Matrix&);
	T operator()(int row, int col) const;

	void setValue(int row, int col, T val);
	void print() const;

	T determinant() const;
	Matrix transpose() const;
	Matrix inverse() const;
	std::vector<Matrix> LU_factorization();
	std::vector<Matrix> QR_factorization(); // with Grand Schmidt method

	const int getNumberOfRows() const 
	{ return numberOfRows; }
	
	const int getNumberOfColumns() const 
	{ return numberOfColumns; }

private:
	T** values;
	int numberOfRows;
	int numberOfColumns;
};

template <class T>
Matrix<T> operator+(const Matrix<T> &A, const Matrix<T> &B)
{
	assert(A.getNumberOfRows() == B.getNumberOfRows());
	assert(A.getNumberOfColumns() == B.getNumberOfColumns());
	Matrix<T> result(A.getNumberOfRows(), A.getNumberOfColumns());
	for (int i = 0; i < A.getNumberOfRows(); ++i)
	{
		for (int j = 0; j < A.getNumberOfColumns(); ++j)
		{
			result.setValue(i, j, A(i,j) + B(i,j));
		}
	}

	return result;
}

template <class T>
Matrix<T> operator-(const Matrix<T> &A, const Matrix<T> &B)
{
	assert(A.getNumberOfRows() == B.getNumberOfRows());
	assert(A.getNumberOfColumns() == B.getNumberOfColumns());
	Matrix<T> result(A.getNumberOfRows(), A.getNumberOfColumns());
	for (int i = 0; i < A.getNumberOfRows(); ++i)
	{
		for (int j = 0; j < A.getNumberOfColumns(); ++j)
		{
			result.setValue(i, j, A(i,j) - B(i,j));
		}
	}
	return result;
}

template <class T>
Matrix<T> operator*(const Matrix<T> &A, const Matrix<T> &B)
{
	assert(B.getNumberOfRows() == A.getNumberOfColumns());
	Matrix<T> result(A.getNumberOfRows(), B.getNumberOfColumns());
	for (int i = 0; i < A.getNumberOfRows(); ++i)
	{
		for (int j = 0; j < B.getNumberOfColumns(); ++j)
		{
			for (int k = 0; k < A.getNumberOfColumns(); ++k)
			{
				result.setValue(i,j, result(i,j) + A(i,k) * B(k,j));
			}
		}
	}
	return result;
}

template <class T>
Matrix<T> operator*(const Matrix<T> &A, const T &constant)
{
	Matrix<T> result(A.getNumberOfRows(), A.getNumberOfColumns());
	for (int i = 0; i < A.getNumberOfRows(); ++i)
	{
		for (int j = 0; j < A.getNumberOfColumns(); ++j)
		{
			result.setValue(i, j, A(i,j) * constant);
		}
	}

	return result;
}

template <class T>
std::vector<T> operator*(const Matrix<T> &A, const std::vector<T> &vec)
{
	assert(vec.size() == A.getNumberOfColumns());
	std::vector<T> result(A.getNumberOfRows());
	for (int i = 0; i < A.getNumberOfRows(); ++i)
	{
		for (int j = 0; j < A.getNumberOfColumns(); ++j)
		{
			result[i] += A(i,j) * vec[j];
		}
	}

	return result;
}

template <class T>
inline T Matrix<T>::operator()(int row, int column) const
{
	assert(row < getNumberOfRows());
	assert(column < getNumberOfColumns());
	return values[row][column];
}

#include "..\src\Matrix.cpp"

#endif