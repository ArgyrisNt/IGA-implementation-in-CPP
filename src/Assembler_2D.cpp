#include <iostream>
#include "..\include\Assembler_2D.h"


std::vector<int> Assembler_2D::computeActiveControlPoints(double g1, double g2)
{
	int span_g1 = spanOfValueInKnotVector(g1, 0);
	int span_g2 = spanOfValueInKnotVector(g2, 1);
	std::vector<int> index, index_x, index_y;
	for (int i = 0; i < getBspline_x().getDegree() + 1; ++i)
	{
		int index_x = span_g1 - getBspline_x().getDegree() + i;
		for (int j = 0; j < getBspline_y().getDegree() + 1; ++j)
		{
			int index_y = span_g2 - getBspline_y().getDegree() + j;
			index.push_back(index_x * (getBspline_y().getNumberOfBasisFunctions()) + index_y);
		}
	}

	return index;
}

Matrix<double> Assembler_2D::Jacobian(const std::vector<int> &indices, const std::vector<double> &dNxNy, const std::vector<double> &NxdNy)
{
	std::vector<double> temp_x(2, 0.0);
	for (int kk = 0; kk < dNxNy.size(); ++kk)
	{
		temp_x[0] += dNxNy[kk] * getControlPoints()[indices[kk]].x;
		temp_x[1] += dNxNy[kk] * getControlPoints()[indices[kk]].y;
	}

	std::vector<double> temp_y(2, 0.0);
	for (int kk = 0; kk < NxdNy.size(); ++kk)
	{
		temp_y[0] += NxdNy[kk] * getControlPoints()[indices[kk]].x;
		temp_y[1] += NxdNy[kk] * getControlPoints()[indices[kk]].y;
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
	int Nx = getDistinctKnots(0).size() - 1;
	int Ny = getDistinctKnots(1).size() - 1;
	for (int ie1 = 0; ie1 < Nx; ++ie1)
	{
		for (int ie2 = 0; ie2 < Ny; ++ie2)
		{
			double left_x = getDistinctKnots(0)[ie1];
			double right_x = getDistinctKnots(0)[ie1 + 1];
			double left_y = getDistinctKnots(1)[ie2];
			double right_y = getDistinctKnots(1)[ie2 + 1];
			Element element(getTrimmingCurve());
			std::vector<Vertex<double>> vertices{Vertex<double>(left_x, left_y), Vertex<double>(left_x, right_y),
												 Vertex<double>(right_x, left_y), Vertex<double>(right_x, right_y)};
			element.setVertices(vertices);
			element.categorise();
			elements.push_back(element);
		}
	}
}

void Assembler_2D::computeStiffnessMatrixAndRightHandSide()
{
	computeTrimmedElements();
	// Assemble stiffness matrix
	int Nx = getDistinctKnots(0).size() - 1;
	int Ny = getDistinctKnots(1).size() - 1;
	Matrix<double> A(getNumberOfBasisFunctions(), getNumberOfBasisFunctions());
	std::vector<double> b(getNumberOfBasisFunctions(), 0.0);
	stiffnessMatrix = A, rightHandSide = b;
	int elementId = -1;
	for (int ie1 = 0; ie1 < Nx; ++ie1)
	{
		int spanX = spanOfValueInKnotVector(getDistinctKnots(0)[ie1], 0);
		for (int ie2 = 0; ie2 < Ny; ++ie2)
		{
			int spanY = spanOfValueInKnotVector(getDistinctKnots(1)[ie2], 1);
			elementId++;
			if (elements[elementId].isTrimmed())
			{
				std::vector<Triangle<double>> triangles = elements[elementId].divideInTriangles();
				for (auto triangle : triangles)
				{
					trimmed_triangles.push_back(triangle);
					double a_x = std::min({triangle.vertex1.x, triangle.vertex2.x, triangle.vertex3.x});
					double b_x = std::max({triangle.vertex1.x, triangle.vertex2.x, triangle.vertex3.x});
					double a_y = std::min({triangle.vertex1.y, triangle.vertex2.y, triangle.vertex3.y});
					double b_y = std::max({triangle.vertex1.y, triangle.vertex2.y, triangle.vertex3.y});
					XGaussPointsAndWeights = GaussPointsAndWeightsTriangle(a_x, b_x);
					YGaussPointsAndWeights = GaussPointsAndWeightsTriangle(a_y, b_y);
					AssistComputeStiffnessMatrixAndRightHandSide(spanX, spanY);
				}
			}
			else
			{
				XGaussPointsAndWeights = GaussPointsAndWeightsQuad(getBspline_x().getDegree() + 3, getDistinctKnots(0)[ie1], getDistinctKnots(0)[ie1 + 1]);
				YGaussPointsAndWeights = GaussPointsAndWeightsQuad(getBspline_y().getDegree() + 3, getDistinctKnots(1)[ie2], getDistinctKnots(1)[ie2 + 1]);
				AssistComputeStiffnessMatrixAndRightHandSide(spanX, spanY);
			}
		}
	}
}

void Assembler_2D::AssistComputeStiffnessMatrixAndRightHandSide(int spanX, int spanY)
{
	for (int il_1 = 0; il_1 < getBspline_x().getDegree() + 1; ++il_1)
	{
		int i1 = spanX - getBspline_x().getDegree() + il_1;
		for (int il_2 = 0; il_2 < getBspline_y().getDegree() + 1; ++il_2)
		{
			int i2 = spanY - getBspline_y().getDegree() + il_2;

			int index = i1 * getBspline_y().getNumberOfBasisFunctions() + i2;
			rightHandSide[index] += computeRightHandSideIntegral(il_1, il_2);

			for (int jl_1 = 0; jl_1 < getBspline_x().getDegree() + 1; ++jl_1)
			{
				int j1 = spanX - getBspline_x().getDegree() + jl_1;
				for (int jl_2 = 0; jl_2 < getBspline_y().getDegree() + 1; ++jl_2)
				{
					int j2 = spanY - getBspline_y().getDegree() + jl_2;

					int index1 = i1 * getBspline_y().getNumberOfBasisFunctions() + i2;
					int index2 = j1 * getBspline_y().getNumberOfBasisFunctions() + j2;
					double newValue = stiffnessMatrix(index1, index2) + computeStiffnessIntegral(il_1, jl_1, il_2, jl_2);
					stiffnessMatrix.setValue(index1, index2, newValue);
				}
			}
		}
	}
}

double Assembler_2D::computeStiffnessIntegral(int il_1, int jl_1, int il_2, int jl_2)
{
	double v = 0.0;
	for (int g1 = 0; g1 < XGaussPointsAndWeights.size(); ++g1)
	{
		for (int g2 = 0; g2 < YGaussPointsAndWeights.size(); ++g2)
		{
			std::pair<std::vector<double>, std::vector<double>> eval_1 = getBspline_x().evaluateAtPoint(XGaussPointsAndWeights[g1].first);
			std::pair<std::vector<double>, std::vector<double>> eval_2 = getBspline_y().evaluateAtPoint(YGaussPointsAndWeights[g2].first);

			std::vector<double> dNxNy = createTensorProduct(eval_1.second, eval_2.first);
			std::vector<double> NxdNy = createTensorProduct(eval_1.first, eval_2.second);
			std::vector<int> indices = computeActiveControlPoints(XGaussPointsAndWeights[g1].first, YGaussPointsAndWeights[g2].first);
			Matrix<double> J = Jacobian(indices, dNxNy, NxdNy);

			double bi_0 = eval_1.first[il_1] * eval_2.first[il_2];
			double bi_x = (J(1, 1) * eval_1.second[il_1] * eval_2.first[il_2] - J(1, 0) * eval_1.first[il_1] * eval_2.second[il_2]);
			double bi_y = (-J(0, 1) * eval_1.second[il_1] * eval_2.first[il_2] + J(0, 0) * eval_1.first[il_1] * eval_2.second[il_2]);
			double bj_x = (J(1, 1) * eval_1.second[jl_1] * eval_2.first[jl_2] - J(1, 0) * eval_1.first[jl_1] * eval_2.second[jl_2]);
			double bj_y = (-J(0, 1) * eval_1.second[jl_1] * eval_2.first[jl_2] + J(0, 0) * eval_1.first[jl_1] * eval_2.second[jl_2]);

			double wvol = XGaussPointsAndWeights[g1].second * YGaussPointsAndWeights[g2].second * fabs(J.determinant());

			v += (bi_x * bj_x + bi_y * bj_y) * wvol / std::pow(J.determinant(),2);
		}
	}

	return v;
}

double Assembler_2D::computeRightHandSideIntegral(int il_1, int il_2)
{
	double v = 0.0;
	for (int g1 = 0; g1 < XGaussPointsAndWeights.size(); ++g1)
	{
		for (int g2 = 0; g2 < YGaussPointsAndWeights.size(); ++g2)
		{
			std::pair<std::vector<double>, std::vector<double>> eval_1 = getBspline_x().evaluateAtPoint(XGaussPointsAndWeights[g1].first);
			std::pair<std::vector<double>, std::vector<double>> eval_2 = getBspline_y().evaluateAtPoint(YGaussPointsAndWeights[g2].first);

			std::vector<double> dNxNy = createTensorProduct(eval_1.second, eval_2.first);
			std::vector<double> NxdNy = createTensorProduct(eval_1.first, eval_2.second);
			std::vector<int> indices = computeActiveControlPoints(XGaussPointsAndWeights[g1].first, YGaussPointsAndWeights[g2].first);
			Matrix<double> J = Jacobian(indices, dNxNy, NxdNy);

			double bi_0 = eval_1.first[il_1] * eval_2.first[il_2];
			double wvol = XGaussPointsAndWeights[g1].second * YGaussPointsAndWeights[g2].second * fabs(J.determinant());

			v += bi_0 * sourceFunction * wvol;
		}
	}

	return v;
}

int Assembler_2D::identifyBoundarySideOfBasisFunction(int i)
{
	int N = getBspline_y().getNumberOfBasisFunctions();

	bool isDirichletOnWestBoundary = ((i >= 0 && i <= N - 1) && (boundaryConditions.west.first == "Dirichlet"));
	bool isDirichletOnSouthBoundary = ((i % N == 0) && (boundaryConditions.south.first == "Dirichlet"));
	bool isDirichletOnNorthBoundary = ((i % N == N - 1) && (boundaryConditions.north.first == "Dirichlet"));
	bool isDirichletOnEastBoundary = ((i >= getNumberOfBasisFunctions() - N) && (boundaryConditions.east.first == "Dirichlet"));

	if (isDirichletOnWestBoundary) return 1;
	else if (isDirichletOnSouthBoundary) return 4;
	else if (isDirichletOnNorthBoundary) return 3;
	else if (isDirichletOnEastBoundary) return 2;

	return 0;
}

void Assembler_2D::computeBoundary()
{
	for (int i = 0; i < getNumberOfBasisFunctions(); ++i)
	{
		int side = identifyBoundarySideOfBasisFunction(i);
		if (side) boundaryBasisFunctions.push_back(std::make_pair(i, side));
	}
}
