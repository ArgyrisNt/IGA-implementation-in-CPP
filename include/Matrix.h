#pragma once

#include <iostream>
#include <vector>

template<class T>
class Matrix
{
public:
	// Constructors
	Matrix(); // Default constructor
	Matrix(int my_rows, int my_cols, T val = 0.0); // Initialize matrix by dimensions and init value
	Matrix(const Matrix&); // Copy constructor

	// Destructor
	~Matrix();

	// Overloaded operators
	Matrix& operator=(Matrix); // Equality between two matrices
	Matrix operator+(Matrix); // Addition of two metrices
	Matrix operator-(Matrix&); // Substraction of two metrices
	Matrix operator*(Matrix&); // Multiplication of two metrices
	Matrix operator*(const T&); // Multiplication of a matrix and a constant
	std::vector<T> operator*(std::vector<T> &); // Multiplication of a matrix and a vector
	T operator()(int row, int col); // Overload operator ()

	//Member functions
	void setValue(int row, int col, T val); // Set a value on a specific position
	T determinant(); // Compute determinant
	void print(); // Print matrix
	Matrix transpose(); // Transpose matrix
	Matrix inverse(); // Inverse matrix

	std::vector<T> Jacobi_iterator(std::vector<T> b, int iters); // Jacobi iterative method to solve a linear system
	std::vector<T> GaussSeidel_iterator(std::vector<T> b, int iters); // Gauss Seidel iterative method to solve a linear system
	std::vector<T> SOR_iterator(std::vector<T> b, int iters, double omega);
	std::vector<T> gradient_iterator(std::vector<T> b, int iters);
	std::vector<T> conjugate_gradient_iterator(std::vector<T> &b, int iters);

	std::vector<Matrix> LU_factorization(); // LU decomposition
	std::vector<T> forward_Euler(std::vector<T> b); // forward Euler method to solve a linear system
	std::vector<T> backward_Euler(std::vector<T> b); // backward Euler method to solve a linear system
	std::vector<Matrix> QR_factorization(); // with Grand Schmidt method

	// Member getter functions
	int getNumberOfRows() { return numberOfRows; }
	int getNumberOfColumns() { return numberOfColumns; }

private:
	// Member variables
	T** values;
	int numberOfRows;
	int numberOfColumns;
};

struct CompareFirst
{
	CompareFirst(int val) : val_(val) {}
	bool operator()(const std::pair<int, char> &elem) const
	{
		return val_ == elem.first;
	}

private:
	int val_;
};
