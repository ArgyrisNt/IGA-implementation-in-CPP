#pragma once

#include <iostream>
#include <vector>

class Matrix
{
public:
	// Constructors
	Matrix();
	Matrix(int my_rows, int my_cols);
	Matrix(int my_rows, int my_cols, double val);
	Matrix(const Matrix&);

	// Overloaded operators
	Matrix& operator=(Matrix);
	Matrix operator+(Matrix&);
	
	Matrix operator-(Matrix&);
	Matrix operator*(Matrix&);
	Matrix operator*(const double&);
	std::vector<double> operator*(std::vector<double>&);
	double operator()(int row, int col);

	//Member functions
	void setValue(int row, int col, double val);
	void calcDet();
	void print();
	Matrix transpose();
	Matrix inverse();
	double norm(std::vector<double>& v);
	std::vector<double> Jacobi_iterator(std::vector<double> b, int iters);
	std::vector<Matrix> LU_factor();
	std::vector<double> forward_Euler(std::vector<double> b);
	std::vector<double> backward_Euler(std::vector<double> b);

	// Destructor
	~Matrix();

	//Member variables
	int rows;
	int cols;
	double** mat;
	double detVal;
};