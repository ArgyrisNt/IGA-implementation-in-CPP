#ifndef H_MATRIX
#define H_MATRIX

#include <iostream>
#include <vector>
#include "..\include\Utilities.h"

template<class T>
class Matrix
{
public:
	Matrix(); 
	Matrix(int my_rows, int my_cols, T val = 0.0);
	Matrix(const Matrix&);

	~Matrix();

	Matrix& operator=(const Matrix&);
	Matrix operator+(Matrix &);
	Matrix operator+(Matrix &&);
	Matrix operator-(Matrix&);
	Matrix operator*(Matrix&);
	Matrix operator*(const T&);
	std::vector<T> operator*(std::vector<T> &);
	T operator()(int row, int col);

	void setValue(int row, int col, T val);

	void print();

	T determinant();
	Matrix transpose();
	Matrix inverse();
	std::vector<Matrix> LU_factorization();
	std::vector<Matrix> QR_factorization(); // with Grand Schmidt method
	std::vector<T> forward_Euler(std::vector<T>& b);
	std::vector<T> backward_Euler(std::vector<T>& b);

	int getNumberOfRows() { return numberOfRows; }
	int getNumberOfColumns() { return numberOfColumns; }

private:
	T** values;
	int numberOfRows;
	int numberOfColumns;
};

#include "..\src\Matrix.cpp"

#endif