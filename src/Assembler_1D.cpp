#include <iostream>
#include <cmath>
#include "..\include\Assembler_1D.h"


Matrix<double> Assembler_1D::Jacobian(double g1, const std::vector<double> &Nx)
{
	int span_g1 = spanOfValueInKnotVector(g1, 0);
	double temp1 = 0.0, temp2 = 0.0;
	Matrix<double> Jacobian(1, 2);
	for (int kk = 0; kk < Nx.size(); ++kk)
	{
		Jacobian.setValue(0, 0, Jacobian(0, 0) + Nx[kk] * getControlPoints()[span_g1 - getBspline_x().getDegree() + kk].x);
		Jacobian.setValue(0, 1, Jacobian(0, 1) + Nx[kk] * getControlPoints()[span_g1 - getBspline_x().getDegree() + kk].y);
	}

	return Jacobian;
}


void Assembler_1D::computeStiffnessMatrixAndRightHandSide()
{
	Matrix<double> A(getNumberOfBasisFunctions(), getNumberOfBasisFunctions());
	std::vector<double> b(getNumberOfBasisFunctions(), 0.0);
	int N = getDistinctKnots(0).size() - 1; // Number of elements
	for (int ie1 = 0; ie1 < N; ++ie1)
	{
		XGaussPointsAndWeights = GaussPointsAndWeightsQuad(getBspline_x().getDegree() + 3, getDistinctKnots(0)[ie1], getDistinctKnots(0)[ie1 + 1]);
		int i_span_1 = spanOfValueInKnotVector(getDistinctKnots(0)[ie1], 0);
		for (int il_1 = 0; il_1 < getBspline_x().getDegree() + 1; ++il_1)
		{
			int i1 = i_span_1 - getBspline_x().getDegree() + il_1;
			b[i1] += computeRightHandSideIntegral(il_1);

			for (int jl_1 = 0; jl_1 < getBspline_x().getDegree() + 1; ++jl_1)
			{
				int j1 = i_span_1 - getBspline_x().getDegree() + jl_1;
				double newValue = A(i1, j1) + computeStiffnessIntegral(il_1, jl_1);
				A.setValue(i1, j1, newValue);
			}
		}
	}
	stiffnessMatrix = A;
	rightHandSide = b;
}

double Assembler_1D::computeStiffnessIntegral(int basisFunction, int trialFunction)
{
	double v = 0.0;
	for (int g1 = 0; g1 < XGaussPointsAndWeights.size(); ++g1)
	{
		std::vector<double> gradVal = getBspline_x().evaluateAtPoint(XGaussPointsAndWeights[g1].first).second;
		Matrix<double> J = Jacobian(XGaussPointsAndWeights[g1].first, gradVal);
		double detJ = sqrt(pow(J(0, 0), 2) + pow(J(0, 1), 2));

		double bi_x = (1.0 / detJ) * gradVal[basisFunction];
		double bj_x = (1.0 / detJ) * gradVal[trialFunction];

		double wvol = XGaussPointsAndWeights[g1].second * fabs(detJ);

		v += (bi_x * bj_x) * wvol;
	}

	return v;
}

double Assembler_1D::computeRightHandSideIntegral(int basisFunction)
{
	double v = 0.0;
	for (int g1 = 0; g1 < XGaussPointsAndWeights.size(); ++g1)
	{
		std::pair<std::vector<double>, std::vector<double>> eval = getBspline_x().evaluateAtPoint(XGaussPointsAndWeights[g1].first);

		Matrix<double> J = Jacobian(XGaussPointsAndWeights[g1].first, eval.second);
		double detJ = sqrt(pow(J(0, 0), 2) + pow(J(0, 1), 2));

		double bi_0 = eval.first[basisFunction];

		double wvol = XGaussPointsAndWeights[g1].second * fabs(detJ);
		v += bi_0 * sourceFunction * wvol;
	}

	return v;
}

void Assembler_1D::computeBoundary()
{
	for (int i = 0; i < getNumberOfBasisFunctions(); ++i)
	{
		bool isDirichletOnWestBoundary = ((i == 0) && (getBoundaryConditions().west.first == "Dirichlet"));
		bool isDirichletOnEastBoundary = ((i == getNumberOfBasisFunctions() - 1) && (getBoundaryConditions().east.first == "Dirichlet"));
		if (isDirichletOnWestBoundary) boundaryBasisFunctions.push_back(std::make_pair(i, 1));
		else if (isDirichletOnEastBoundary) boundaryBasisFunctions.push_back(std::make_pair(i, 2));
	}
}
