#include <iostream>
#include "..\include\DiffusionAssembler_1D.h"

Matrix<double> &DiffusionAssembler_1D::getMassMatrix()
{
	return massMatrix;
}



void DiffusionAssembler_1D::assemble()
{
	computeStiffnessMatrixAndRightHandSide();
	computeMassMatrix();
	computeBoundary();
	systemMatrix = massMatrix + stiffnessMatrix * (coefficient * Timestep);
	for (int i = 0; i < rightHandSide.size(); ++i)
	{
		rightHandSide[i] = rightHandSide[i] * Timestep;
	}
}

void DiffusionAssembler_1D::computeMassMatrix()
{
	// Assemble mass matrix
	Matrix<double> M(getNumberOfBasisFunctions(), getNumberOfBasisFunctions());
	int N = getDistinctKnots(0).size() - 1; // Number of elements
	for (int ie1 = 0; ie1 < N; ++ie1)
	{
		int i_span_1 = spanOfValueInKnotVector(getDistinctKnots(0)[ie1], 0);
		for (int il_1 = 0; il_1 < getBspline_x().getDegree() + 1; ++il_1)
		{
			int i1 = i_span_1 - getBspline_x().getDegree() + il_1;
			for (int jl_1 = 0; jl_1 < getBspline_x().getDegree() + 1; ++jl_1)
			{
				int j1 = i_span_1 - getBspline_x().getDegree() + jl_1;

				double newValue = M(i1, j1) + computeMassIntegral(ie1, il_1, jl_1);
				M.setValue(i1, j1, newValue);
			}
		}
	}
	massMatrix = M;
}

double DiffusionAssembler_1D::computeMassIntegral(int element, int basisFunction, int trialFunction)
{
	double v = 0.0;
	std::vector<std::pair<double, double>> gauss = GaussPointsAndWeightsQuad(getBspline_x().getDegree() + 3, getDistinctKnots(0)[element], getDistinctKnots(0)[element + 1]);
	for (int g1 = 0; g1 < gauss.size(); ++g1)
	{
		std::pair<std::vector<double>, std::vector<double>> eval = getBspline_x().evaluateAtPoint(gauss[g1].first);
		Matrix<double> J = Jacobian(gauss[g1].first, eval.second);
		double detJ = sqrt(pow(J(0, 0), 2) + pow(J(0, 1), 2));

		double bi_0 = eval.first[basisFunction];
		double bj_0 = eval.first[trialFunction];

		double wvol = gauss[g1].second * fabs(detJ);

		v += (bi_0 * bj_0) * wvol;
	}
	return v;
}

std::vector<double> DiffusionAssembler_1D::nextStep(const std::vector<double> &sol)
{
    std::vector<double> b(rightHandSide.size(), 0.0);
	std::vector<double> temp = massMatrix * sol;
	for (int i = 0; i < rightHandSide.size(); ++i)
	{
		b[i] = rightHandSide[i] + temp[i];
	}

	return b;
}

std::vector<double> DiffusionAssembler_1D::applyInitialCondition(double (*func)(double))
{
	std::vector<double> sol;
	std::vector<double> linspaced;
	double start = getDistinctKnots(0)[0];
	double end = getDistinctKnots(0)[getDistinctKnots(0).size() - 1];
	double delta = (end - start) / (getNumberOfBasisFunctions() - 1);
	for (int i = 0; i < getNumberOfBasisFunctions() - 1; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end);
	for (int i = 0; i < linspaced.size(); ++i)
	{
		sol.push_back(func(linspaced[i]));
	}
	return sol;
}
