#include <iostream>
#include "..\include\Assembler_2D.h"

Bspline &Assembler_2D::getBspline_y()
{
	return bsplineEntity.getBspline(1);
}

Bspline &Assembler_2D::getBspline_x()
{
	return bsplineEntity.getBspline(0);
}

std::vector<std::vector<double>> &Assembler_2D::getControlPoints()
{
	return bsplineEntity.getControlPoints();
}

const int Assembler_2D::getNumberOfBasisFunctions()
{
	return getBspline_x().getNumberOfBasisFunctions() * getBspline_y().getNumberOfBasisFunctions();
}

TrimmingCurve &Assembler_2D::getTrimmingCurve()
{ 
	return bsplineEntity.trimmingCurve;
}



std::vector<int> Assembler_2D::computeActiveControlPoints(double g1, double g2)
{
	int span_g1 = spanOfValueInKnotVector(g1, 0);
	int span_g2 = spanOfValueInKnotVector(g2, 1);
	std::vector<int> index, index_x, index_y;
	for (int kk = 0; kk < getBspline_x().getDegree() + 1; kk++)
	{
		index_x.push_back(span_g1 - getBspline_x().getDegree() + kk);
	}
	for (int kk = 0; kk < getBspline_y().getDegree() + 1; kk++)
	{
		index_y.push_back(span_g2 - getBspline_y().getDegree() + kk);
	}
	for (auto indx : index_x)
	{
		for (auto indy : index_y)
		{
			index.push_back(indx * (getBspline_y().getNumberOfBasisFunctions()) + indy);
		}
	}
	return index;
}

Matrix<double> Assembler_2D::Jacobian(std::vector<int>& indices, std::vector<double> &dNxNy, std::vector<double> &NxdNy)
{
	std::vector<double> temp_x(2, 0.0);
	for (int kk = 0; kk < dNxNy.size(); kk++)
	{
		temp_x[0] += dNxNy[kk] * getControlPoints()[indices[kk]][0];
		temp_x[1] += dNxNy[kk] * getControlPoints()[indices[kk]][1];
	}

	std::vector<double> temp_y(2, 0.0);
	for (int kk = 0; kk < NxdNy.size(); kk++)
	{
		temp_y[0] += NxdNy[kk] * getControlPoints()[indices[kk]][0];
		temp_y[1] += NxdNy[kk] * getControlPoints()[indices[kk]][1];
	}

	Matrix<double> Jacobian(2, 2);
	Jacobian.setValue(0, 0, temp_x[0]);
	Jacobian.setValue(0, 1, temp_y[0]);
	Jacobian.setValue(1, 0, temp_x[1]);
	Jacobian.setValue(1, 1, temp_y[1]);

	return Jacobian;
}



void Assembler_2D::assemble()
{
	computeTrimmedElements();
	computeStiffnessMatrixAndRighHandSide();
	computeBoundary();
	systemMatrix = stiffnessMatrix;
}

void Assembler_2D::computeTrimmedElements()
{
	int Nx = getDistinctKnots(0).size() - 1;
	int Ny = getDistinctKnots(1).size() - 1;
	for (int ie1 = 0; ie1 < Nx; ie1++)
	{
		for (int ie2 = 0; ie2 < Ny; ie2++)
		{
			double left_x = getDistinctKnots(0)[ie1];
			double right_x = getDistinctKnots(0)[ie1 + 1];
			double left_y = getDistinctKnots(1)[ie2];
			double right_y = getDistinctKnots(1)[ie2 + 1];
			Element element(false, getTrimmingCurve());
			std::vector<Vertex<double>> vertices{Vertex<double>(left_x, left_y), Vertex<double>(left_x, right_y),
												 Vertex<double>(right_x, left_y), Vertex<double>(right_x, right_y)};
			element.setVertices(vertices);
			element.categorise();
			elements.push_back(element);
		}
	}
}

void Assembler_2D::computeStiffnessMatrixAndRighHandSide()
{
	// Assemble stiffness matrix
	int Nx = getDistinctKnots(0).size() - 1;
	int Ny = getDistinctKnots(1).size() - 1;
	Matrix<double> A(getNumberOfBasisFunctions(), getNumberOfBasisFunctions());
	std::vector<double> b(getNumberOfBasisFunctions(), 0.0);
	stiffnessMatrix = A, rightHandSide = b;
	int elementId = -1;
	for (int ie1 = 0; ie1 < Nx; ie1++)
	{
		for (int ie2 = 0; ie2 < Ny; ie2++)
		{
			elementId++;
			if (elements[elementId].isTrimmed)
			{
				std::vector<Triangle<double>> triangles = elements[elementId].divideInTriangles();
				for (auto triangle : triangles)
				{
					trimmed_triangles.push_back(triangle);
					computeTriangleStiffnessMatrixAndRightHandSide(triangle, ie1, ie2);
				}
			}
			else
			{
				computeQuadStiffnessMatrixAndRightHandSide(ie1, ie2);
			}
		}
	}
}

void Assembler_2D::computeQuadStiffnessMatrixAndRightHandSide(int elementX, int elementY)
{
	int spanX = spanOfValueInKnotVector(getDistinctKnots(0)[elementX], 0);
	int spanY = spanOfValueInKnotVector(getDistinctKnots(1)[elementY], 1);
	XGaussPointsAndWeights = GaussPointsAndWeightsQuad(getBspline_x().getDegree() + 3, getDistinctKnots(0)[elementX], getDistinctKnots(0)[elementX + 1]);
	YGaussPointsAndWeights = GaussPointsAndWeightsQuad(getBspline_y().getDegree() + 3, getDistinctKnots(1)[elementY], getDistinctKnots(1)[elementY + 1]);
	for (int il_1 = 0; il_1 < getBspline_x().getDegree() + 1; il_1++)
	{
		int i1 = spanX - getBspline_x().getDegree() + il_1;
		for (int il_2 = 0; il_2 < getBspline_y().getDegree() + 1; il_2++)
		{
			int i2 = spanY - getBspline_y().getDegree() + il_2;

			int index = i1 * getBspline_y().getNumberOfBasisFunctions() + i2;
			rightHandSide[index] += computeRightHandSideIntegral(il_1, il_2);

			for (int jl_1 = 0; jl_1 < getBspline_x().getDegree() + 1; jl_1++)
			{
				int j1 = spanX - getBspline_x().getDegree() + jl_1;
				for (int jl_2 = 0; jl_2 < getBspline_y().getDegree() + 1; jl_2++)
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

void Assembler_2D::computeTriangleStiffnessMatrixAndRightHandSide(Triangle<double> &triangle, int elementX, int elementY)
{
	int spanX = spanOfValueInKnotVector(getDistinctKnots(0)[elementX], 0);
	int spanY = spanOfValueInKnotVector(getDistinctKnots(1)[elementY], 1);
	trimmed_triangles.push_back(triangle);
	double a_x = std::min({triangle.vertex1.x, triangle.vertex2.x, triangle.vertex3.x});
	double b_x = std::max({triangle.vertex1.x, triangle.vertex2.x, triangle.vertex3.x});
	double a_y = std::min({triangle.vertex1.y, triangle.vertex2.y, triangle.vertex3.y});
	double b_y = std::max({triangle.vertex1.y, triangle.vertex2.y, triangle.vertex3.y});
	XGaussPointsAndWeights = GaussPointsAndWeightsTriangle(a_x, b_x);
	YGaussPointsAndWeights = GaussPointsAndWeightsTriangle(a_y, b_y);
	for (int il_1 = 0; il_1 < getBspline_x().getDegree() + 1; il_1++)
	{
		int i1 = spanX - getBspline_x().getDegree() + il_1;
		for (int il_2 = 0; il_2 < getBspline_y().getDegree() + 1; il_2++)
		{
			int i2 = spanY - getBspline_y().getDegree() + il_2;

			int index = i1 * getBspline_y().getNumberOfBasisFunctions() + i2;
			rightHandSide[index] += computeRightHandSideIntegral(il_1, il_2);

			for (int jl_1 = 0; jl_1 < getBspline_x().getDegree() + 1; jl_1++)
			{
				int j1 = spanX - getBspline_x().getDegree() + jl_1;
				for (int jl_2 = 0; jl_2 < getBspline_y().getDegree() + 1; jl_2++)
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
	for (int g1 = 0; g1 < XGaussPointsAndWeights.size(); g1++)
	{
		for (int g2 = 0; g2 < YGaussPointsAndWeights.size(); g2++)
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
	for (int g1 = 0; g1 < XGaussPointsAndWeights.size(); g1++)
	{
		for (int g2 = 0; g2 < YGaussPointsAndWeights.size(); g2++)
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

	bool isDirichletOnWestBoundary = ((i >= 0 && i <= N - 1) && (boundaryConditions->west.first == "Dirichlet"));
	bool isDirichletOnSouthBoundary = ((i % N == 0) && (boundaryConditions->south.first == "Dirichlet"));
	bool isDirichletOnNorthBoundary = ((i % N == N - 1) && (boundaryConditions->north.first == "Dirichlet"));
	bool isDirichletOnEastBoundary = ((i >= getNumberOfBasisFunctions() - N) && (boundaryConditions->east.first == "Dirichlet"));

	if (isDirichletOnWestBoundary && IDHasNotAlreadyBeingMarked(boundaryBasisFunctions, i)) return 1;
	else if (isDirichletOnSouthBoundary && IDHasNotAlreadyBeingMarked(boundaryBasisFunctions, i)) return 4;
	else if (isDirichletOnNorthBoundary && IDHasNotAlreadyBeingMarked(boundaryBasisFunctions, i)) return 3;
	else if (isDirichletOnEastBoundary && IDHasNotAlreadyBeingMarked(boundaryBasisFunctions, i)) return 2;

	return 0;
}

void Assembler_2D::computeBoundary()
{
	for (int i = 0; i < getNumberOfBasisFunctions(); i++)
	{
		int side = identifyBoundarySideOfBasisFunction(i);
		if (side) boundaryBasisFunctions.push_back(std::make_pair(i, side));
	}
}



void Assembler_2D::writeParameterSpaceToFile(std::string filename)
{
	std::ofstream my_file(filename);
	my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
	my_file << "zone t= " << "\"1\"" << ",i=" << getDistinctKnots(0).size() << ",j=" << getDistinctKnots(1).size() << "\n";
	for (int j = 0; j < getDistinctKnots(1).size(); j++)
	{
		for (int i = 0; i < getDistinctKnots(0).size(); i++)
		{
			my_file << getDistinctKnots(0)[i] << " " << getDistinctKnots(1)[j] << "\n";
		}
	}
	my_file.close();
}

void Assembler_2D::writeTrimmedTrianglesToFile(std::string filename)
{
	std::ofstream my_file(filename);
	my_file << "# " << trimmed_triangles.size() << "\n";
	int cnt = 0;
	for (auto triangle : trimmed_triangles)
	{
		cnt += 2;
		my_file << "v " << triangle.vertex1.x << " " << triangle.vertex1.y << "\n";
		my_file << "v " << triangle.vertex2.x << " " << triangle.vertex2.y << "\n";
		my_file << "l " << cnt - 1 << " " << cnt << "\n";
		cnt += 2;
		my_file << "v " << triangle.vertex2.x << " " << triangle.vertex2.y << "\n";
		my_file << "v " << triangle.vertex3.x << " " << triangle.vertex3.y << "\n";
		my_file << "l " << cnt - 1 << " " << cnt << "\n";
		cnt += 2;
		my_file << "v " << triangle.vertex3.x << " " << triangle.vertex3.y << "\n";
		my_file << "v " << triangle.vertex1.x << " " << triangle.vertex1.y << "\n";
		my_file << "l " << cnt - 1 << " " << cnt << "\n";
	}
	my_file.close();
}