#include <iostream>
#include "..\include\Assembler_2D.h"

bool Assembler_2D::basisFunctionHasNotAlreadyMarked(int id)
{
	auto it = std::find_if(boundaryBasisFunctions.begin(), boundaryBasisFunctions.end(), CompareFirst(id));
	return it == boundaryBasisFunctions.end();
}

void Assembler_2D::computeBoundary()
{
	for (int i = 0; i < getNumberOfBasisFunctions(); i++)
	{
		bool isDirichletOnWestBoundary = ((i >= 0 && i <= bspline_y->getNumberOfBasisFunctions() - 1) && (boundaryConditions->getWestType() == "Dirichlet"));
		bool isDirichletOnSouthBoundary = ((i % bspline_y->getNumberOfBasisFunctions() == 0) && (boundaryConditions->getSouthType() == "Dirichlet"));
		bool isDirichletOnNorthBoundary = ((i % bspline_y->getNumberOfBasisFunctions() == bspline_y->getNumberOfBasisFunctions() - 1) && (boundaryConditions->getNorthType() == "Dirichlet"));
		bool isDirichletOnEastBoundary = ((i >= getNumberOfBasisFunctions() - bspline_y->getNumberOfBasisFunctions()) && (boundaryConditions->getEastType() == "Dirichlet"));

		if (isDirichletOnWestBoundary && basisFunctionHasNotAlreadyMarked(i)) boundaryBasisFunctions.push_back(std::make_pair(i, 1));
		else if (isDirichletOnSouthBoundary && basisFunctionHasNotAlreadyMarked(i)) boundaryBasisFunctions.push_back(std::make_pair(i, 4));
		else if (isDirichletOnNorthBoundary && basisFunctionHasNotAlreadyMarked(i)) boundaryBasisFunctions.push_back(std::make_pair(i, 3));
		else if (isDirichletOnEastBoundary && basisFunctionHasNotAlreadyMarked(i)) boundaryBasisFunctions.push_back(std::make_pair(i, 2));
	}
}

std::vector<double> Assembler_2D::createTensorProduct(std::vector<double>& vec1, std::vector<double>& vec2)
{
	std::vector<double> tensor_product;
	for (int u = 0; u < vec1.size(); u++)
	{
		for (int v = 0; v < vec2.size(); v++)
		{
			tensor_product.push_back({ vec1[u] * vec2[v] });
		}
	}

	return tensor_product;
}

std::vector<int> Assembler_2D::computeActiveControlPoints(double g1, double g2)
{
	int span_g1 = bspline_x->getKnotvector().findSpanOfValue(g1);
	int span_g2 = bspline_y->getKnotvector().findSpanOfValue(g2);
	std::vector<int> index, index_x, index_y;
	for (int kk = 0; kk < bspline_x->getDegree() + 1; kk++)
	{
		index_x.push_back(span_g1 - bspline_x->getDegree() + kk);
	}
	for (int kk = 0; kk < bspline_y->getDegree() + 1; kk++)
	{
		index_y.push_back(span_g2 - bspline_y->getDegree() + kk);
	}
	for (auto indx : index_x)
	{
		for (auto indy : index_y)
		{
			index.push_back(indx * (bspline_y->getNumberOfBasisFunctions()) + indy);
		}
	}
	return index;
}

Matrix<double> Assembler_2D::Jacobian(double g1, double g2, std::pair<std::vector<double>, std::vector<double>> &eval_1, std::pair<std::vector<double>, std::vector<double>>& eval_2)
{	
	std::vector<double> shp_fnc_dx = createTensorProduct(eval_1.second, eval_2.first);
	std::vector<double> shp_fnc_dy = createTensorProduct(eval_1.first, eval_2.second);

	std::vector<int> index = computeActiveControlPoints(g1, g2);
	
	std::vector<double>	temp_x(2, 0.0);
	for (int kk = 0; kk < shp_fnc_dx.size(); kk++)
	{
		temp_x[0] += shp_fnc_dx[kk] * controlPoints[index[kk]][0];
		temp_x[1] += shp_fnc_dx[kk] * controlPoints[index[kk]][1];
	}

	std::vector<double>	temp_y(2, 0.0);
	for (int kk = 0; kk < shp_fnc_dy.size(); kk++)
	{
		temp_y[0] += shp_fnc_dy[kk] * controlPoints[index[kk]][0];
		temp_y[1] += shp_fnc_dy[kk] * controlPoints[index[kk]][1];
	}

	Matrix<double> Jacobian(2, 2);
	Jacobian.setValue(0, 0, temp_x[0]);
	Jacobian.setValue(0, 1, temp_y[0]);
	Jacobian.setValue(1, 0, temp_x[1]);
	Jacobian.setValue(1, 1, temp_y[1]);

	return Jacobian;
}

void Assembler_2D::computeTrimmedElements()
{
	int Nx = bspline_x->getKnotvector().getDistinctKnots().size() - 1;
	int Ny = bspline_y->getKnotvector().getDistinctKnots().size() - 1;
	for (int ie1 = 0; ie1 < Nx; ie1++)
	{
		for (int ie2 = 0; ie2 < Ny; ie2++)
		{
			double left_x = bspline_x->getKnotvector().getDistinctKnots()[ie1];
			double right_x = bspline_x->getKnotvector().getDistinctKnots()[ie1 + 1];
			double left_y = bspline_y->getKnotvector().getDistinctKnots()[ie2];
			double right_y = bspline_y->getKnotvector().getDistinctKnots()[ie2 + 1];
			Element element(false, trimmingCurve);
			element.vertices = {std::make_pair(left_x, left_y), std::make_pair(left_x, right_y), std::make_pair(right_x, left_y), std::make_pair(right_x, right_y)};
			element.categorise();
			elements.push_back(element);
		}
	}
}

void Assembler_2D::computeStiffnessMatrixAndRighHandSide()
{
	// Assemble stiffness matrix
	int Nx = bspline_x->getKnotvector().getDistinctKnots().size() - 1; // Number of elements on x-direction
	int Ny = bspline_y->getKnotvector().getDistinctKnots().size() - 1; // NUmber of elements on y-direction
	Matrix<double> A(getNumberOfBasisFunctions(), getNumberOfBasisFunctions());
	std::vector<double> b(getNumberOfBasisFunctions(), 0.0);
	int elementId = -1;
	for (int ie1 = 0; ie1 < Nx; ie1++)
	{
		for (int ie2 = 0; ie2 < Ny; ie2++)
		{
			elementId++;
			if (elements[elementId].isTrimmed)
			{
				std::vector<std::vector<std::pair<double, double>>> triangles = elements[elementId].divideInTriangles();
				for (auto triangle : triangles)
				{
					trimmed_triangles.push_back(triangle);
					computeTriangleStiffnessMatrixAndRightHandSide(triangle, ie1, ie2, A, b);
				}
			}
			else
			{
				computeQuadStiffnessMatrixAndRightHandSide(ie1, ie2, A, b);
			}
		}
	}
	stiffnessMatrix = A;
	rightHandSide = b;
}

void Assembler_2D::computeQuadStiffnessMatrixAndRightHandSide(int elementX, int elementY, Matrix<double> &A, std::vector<double> &b)
{
	int spanX = bspline_x->getKnotvector().findSpanOfValue(bspline_x->getKnotvector().getDistinctKnots()[elementX]);
	int spanY = bspline_y->getKnotvector().findSpanOfValue(bspline_y->getKnotvector().getDistinctKnots()[elementY]);
	std::pair<std::vector<double>, std::vector<double>> gauss_x = bspline_x->GaussPointsAndWeights(bspline_x->getDegree() + 3, bspline_x->getKnotvector().getDistinctKnots()[elementX], bspline_x->getKnotvector().getDistinctKnots()[elementX + 1]);
	std::pair<std::vector<double>, std::vector<double>> gauss_y = bspline_y->GaussPointsAndWeights(bspline_y->getDegree() + 3, bspline_y->getKnotvector().getDistinctKnots()[elementY], bspline_y->getKnotvector().getDistinctKnots()[elementY + 1]);
	for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
	{
		int i1 = spanX - bspline_x->getDegree() + il_1;
		for (int il_2 = 0; il_2 < bspline_y->getDegree() + 1; il_2++)
		{
			int i2 = spanY - bspline_y->getDegree() + il_2;

			int index = i1 * bspline_y->getNumberOfBasisFunctions() + i2;
			b[index] += computeRightHandSideIntegral(gauss_x, il_1, gauss_y, il_2);

			for (int jl_1 = 0; jl_1 < bspline_x->getDegree() + 1; jl_1++)
			{
				int j1 = spanX - bspline_x->getDegree() + jl_1;
				for (int jl_2 = 0; jl_2 < bspline_y->getDegree() + 1; jl_2++)
				{
					int j2 = spanY - bspline_y->getDegree() + jl_2;

					int index1 = i1 * bspline_y->getNumberOfBasisFunctions() + i2;
					int index2 = j1 * bspline_y->getNumberOfBasisFunctions() + j2;
					double newValue = A(index1, index2) + computeStiffnessIntegral(gauss_x, il_1, jl_1, gauss_y, il_2, jl_2);
					A.setValue(index1, index2, newValue);
				}
			}
		}
	}
}

void Assembler_2D::computeTriangleStiffnessMatrixAndRightHandSide(std::vector<std::pair<double, double>> &triangle, int elementX, int elementY, Matrix<double> &A, std::vector<double> &b)
{
	int spanX = bspline_x->getKnotvector().findSpanOfValue(bspline_x->getKnotvector().getDistinctKnots()[elementX]);
	int spanY = bspline_y->getKnotvector().findSpanOfValue(bspline_y->getKnotvector().getDistinctKnots()[elementY]);
	trimmed_triangles.push_back(triangle);
	double a_x = std::min({triangle[0].first, triangle[1].first, triangle[2].first});
	double b_x = std::max({triangle[0].first, triangle[1].first, triangle[2].first});
	double a_y = std::min({triangle[0].second, triangle[1].second, triangle[2].second});
	double b_y = std::max({triangle[0].second, triangle[1].second, triangle[2].second});
	std::pair<std::vector<double>, std::vector<double>> gauss_x = GaussPointsAndWeightsTria(a_x, b_x);
	std::pair<std::vector<double>, std::vector<double>> gauss_y = GaussPointsAndWeightsTria(a_y, b_y);
	for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
	{
		int i1 = spanX - bspline_x->getDegree() + il_1;
		for (int il_2 = 0; il_2 < bspline_y->getDegree() + 1; il_2++)
		{
			int i2 = spanY - bspline_y->getDegree() + il_2;

			int index = i1 * bspline_y->getNumberOfBasisFunctions() + i2;
			b[index] += computeRightHandSideIntegral(gauss_x, il_1, gauss_y, il_2);

			for (int jl_1 = 0; jl_1 < bspline_x->getDegree() + 1; jl_1++)
			{
				int j1 = spanX - bspline_x->getDegree() + jl_1;
				for (int jl_2 = 0; jl_2 < bspline_y->getDegree() + 1; jl_2++)
				{
					int j2 = spanY - bspline_y->getDegree() + jl_2;

					int index1 = i1 * bspline_y->getNumberOfBasisFunctions() + i2;
					int index2 = j1 * bspline_y->getNumberOfBasisFunctions() + j2;
					double newValue = A(index1, index2) + computeStiffnessIntegral(gauss_x, il_1, jl_1, gauss_y, il_2, jl_2);
					A.setValue(index1, index2, newValue);
				}
			}
		}
	}
}

double Assembler_2D::computeStiffnessIntegral(std::pair<std::vector<double>, std::vector<double>> &gauss_x, int il_1, int jl_1, std::pair<std::vector<double>, std::vector<double>> &gauss_y, int il_2, int jl_2)
{
	double v = 0.0;
	for (int g1 = 0; g1 < gauss_x.first.size(); g1++)
	{
		for (int g2 = 0; g2 < gauss_y.first.size(); g2++)
		{
			std::pair<std::vector<double>, std::vector<double>> eval_1 = bspline_x->evaluateAtPoint(gauss_x.first[g1]);
			std::pair<std::vector<double>, std::vector<double>> eval_2 = bspline_y->evaluateAtPoint(gauss_y.first[g2]);

			Matrix<double> J = Jacobian(gauss_x.first[g1], gauss_y.first[g2], eval_1, eval_2);

			double bi_0 = eval_1.first[il_1] * eval_2.first[il_2];
			double bi_x = (1.0 / J.determinant()) * (J(1, 1) * eval_1.second[il_1] * eval_2.first[il_2] - J(1, 0) * eval_1.first[il_1] * eval_2.second[il_2]);
			double bi_y = (1.0 / J.determinant()) * (-J(0, 1) * eval_1.second[il_1] * eval_2.first[il_2] + J(0, 0) * eval_1.first[il_1] * eval_2.second[il_2]);
			double bj_x = (1.0 / J.determinant()) * (J(1, 1) * eval_1.second[jl_1] * eval_2.first[jl_2] - J(1, 0) * eval_1.first[jl_1] * eval_2.second[jl_2]);
			double bj_y = (1.0 / J.determinant()) * (-J(0, 1) * eval_1.second[jl_1] * eval_2.first[jl_2] + J(0, 0) * eval_1.first[jl_1] * eval_2.second[jl_2]);

			double wvol = gauss_x.second[g1] * gauss_y.second[g2] * fabs(J.determinant());

			v += (bi_x * bj_x + bi_y * bj_y) * wvol;
		}
	}

	return v;
}

double Assembler_2D::computeRightHandSideIntegral(std::pair<std::vector<double>, std::vector<double>> &gauss_x, int il_1, std::pair<std::vector<double>, std::vector<double>> &gauss_y, int il_2)
{
	double v = 0.0;
	for (int g1 = 0; g1 < gauss_x.first.size(); g1++)
	{
		for (int g2 = 0; g2 < gauss_y.first.size(); g2++)
		{
			std::pair<std::vector<double>, std::vector<double>> eval_1 = bspline_x->evaluateAtPoint(gauss_x.first[g1]);
			std::pair<std::vector<double>, std::vector<double>> eval_2 = bspline_y->evaluateAtPoint(gauss_y.first[g2]);

			Matrix<double> J = Jacobian(gauss_x.first[g1], gauss_y.first[g2], eval_1, eval_2);

			double bi_0 = eval_1.first[il_1] * eval_2.first[il_2];
			double wvol = gauss_x.second[g1] * gauss_y.second[g2] * fabs(J.determinant());

			v += bi_0 * sourceFunction * wvol;
		}
	}

	return v;
}

std::pair<std::vector<double>, std::vector<double>> Assembler_2D::GaussPointsAndWeightsTria(double a, double b)
{
	std::vector<double> GS_pts_temp{1.0/6.0, 2.0/3.0};
	std::vector<double> GS_wgts{1.0/3.0, 1.0/3.0};
	// Convert to interval [a,b]
	std::vector<double> GS_pts;
	for (int i = 0; i < GS_pts_temp.size(); i++)
	{
		GS_pts.push_back(a * (1 - GS_pts_temp[i]) + b * GS_pts_temp[i]);
	}

	return std::make_pair(GS_pts, GS_wgts);
}
