#include <iostream>
#include "..\include\DiffusionAssembler_1D.h"

double DiffusionAssembler_1D::computeMassIntegral(int element, int basisFunction, int trialFunction)
{
	double v = 0.0;
	std::vector<std::pair<double,double>> gauss = GaussPointsAndWeights(bspline_x->getDegree() + 3, XdistinctKnots[element], XdistinctKnots[element + 1]);
	for (int g1 = 0; g1 < gauss.size(); g1++)
	{
		std::pair<std::vector<double>, std::vector<double>> eval = bspline_x->evaluateAtPoint(gauss[g1].first);
		Matrix<double> J = Jacobian(gauss[g1].first, eval.second);
		double detJ = sqrt(pow(J(0, 0), 2) + pow(J(0, 1), 2));

		double bi_0 = eval.first[basisFunction];
		double bj_0 = eval.first[trialFunction];

		double wvol = gauss[g1].second * fabs(detJ);

		v += (bi_0 * bj_0) * wvol;
	}
	return v;
}

void DiffusionAssembler_1D::computeMassMatrix()
{
    // Assemble mass matrix
	Matrix<double> M(getNumberOfBasisFunctions(), getNumberOfBasisFunctions());
	int N = XdistinctKnots.size() - 1; // Number of elements
	for (int ie1 = 0; ie1 < N; ie1++)
	{
		int i_span_1 = XspanOfValueInKnotVector(XdistinctKnots[ie1]);
		for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
		{
			int i1 = i_span_1 - bspline_x->getDegree() + il_1;
			for (int jl_1 = 0; jl_1 < bspline_x->getDegree() + 1; jl_1++)
			{
				int j1 = i_span_1 - bspline_x->getDegree() + jl_1;

				double newValue = M(i1, j1) + computeMassIntegral(ie1, il_1, jl_1);
				M.setValue(i1, j1, newValue);
			}
		}
	}
	massMatrix = M;
}

std::vector<double> DiffusionAssembler_1D::nextStep(std::vector<double> sol)
{
    std::vector<double> b(rightHandSide.size(), 0.0);
	std::vector<double> temp = massMatrix * sol;
	for (int i = 0; i < rightHandSide.size(); i++)
	{
		b[i] = rightHandSide[i] + temp[i];
	}

	return b;
}

std::vector<double> DiffusionAssembler_1D::applyInitialCondition(double (*func)(double))
{
	std::vector<double> sol;
	std::vector<double> linspaced;
	double start = XdistinctKnots[0];
	double end = XdistinctKnots[XdistinctKnots.size() - 1];
	double delta = (end - start) / (getNumberOfBasisFunctions() - 1);
	for (int i = 0; i < getNumberOfBasisFunctions() - 1; i++)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end);
	for (int i = 0; i < linspaced.size(); i++)
	{
		sol.push_back(func(linspaced[i]));
	}
	return sol;
}

void DiffusionAssembler_1D::applyBoundaryEllimination()
{
	int newDimension = systemMatrix.getNumberOfRows() - boundaryBasisFunctions.size();
	Matrix<double> newSystemMatrix(newDimension, newDimension);
	std::vector<double> newRightHandSide;
	int i = 0;
	for (int ii = 0; ii < stiffnessMatrix.getNumberOfRows(); ii++)
	{
		auto it = std::find_if(boundaryBasisFunctions.begin(), boundaryBasisFunctions.end(), CompareFirst(ii));
		if (it != boundaryBasisFunctions.end()) continue;
		int j = 0;
		for (int jj = 0; jj < stiffnessMatrix.getNumberOfRows(); jj++)
		{
			it = std::find_if(boundaryBasisFunctions.begin(), boundaryBasisFunctions.end(), CompareFirst(jj));
			if (it != boundaryBasisFunctions.end()) continue;
			newSystemMatrix.setValue(i, j, systemMatrix(ii, jj));
			j++;
		}
		newRightHandSide.push_back(rightHandSide[ii]);
		i++;
	}

	systemMatrix = newSystemMatrix;
	rightHandSide = newRightHandSide;
}

void DiffusionAssembler_1D::applyBoundaryMultipliers()
{
	int new_dim = systemMatrix.getNumberOfRows() + boundaryBasisFunctions.size();
	Matrix<double> newSystemMatrix(new_dim, new_dim, 0.0);
	std::vector<double> newRightHandSide;
	int cnt2 = 0;
	for (auto it = boundaryBasisFunctions.begin(); it != boundaryBasisFunctions.end(); it++)
	{
		newSystemMatrix.setValue(cnt2, (*it).first + boundaryBasisFunctions.size(), 1.0);
		newSystemMatrix.setValue((*it).first + boundaryBasisFunctions.size(), cnt2, 1.0);
		newRightHandSide.push_back(addBoundaryValueToRhs(cnt2));
		cnt2++;
	}

	int ii = 0;
	for (int i = cnt2; i < newSystemMatrix.getNumberOfRows(); i++)
	{
		int jj = 0;
		for (int j = cnt2; j < newSystemMatrix.getNumberOfColumns(); j++)
		{
			newSystemMatrix.setValue(i, j, systemMatrix(ii, jj));
			jj++;
		}
		newRightHandSide.push_back(rightHandSide[ii]);
		ii++;
	}

	systemMatrix = newSystemMatrix;
	rightHandSide = newRightHandSide;
}

void DiffusionAssembler_1D::enforceBoundaryConditions(std::string& mode)
{
    boundaryMode = mode;
    if (mode == "Ellimination")
	{	
		applyBoundaryEllimination();
	}
	else if (mode == "Multipliers")
	{
		applyBoundaryMultipliers();
	}
	else
	{	
		std::cout << "Invalid method for enforcing boundary conditions." << std::endl;
		throw std::invalid_argument("Invalid method");
	}
}