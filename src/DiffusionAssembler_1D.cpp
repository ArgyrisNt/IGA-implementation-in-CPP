#include <iostream>
#include "..\include\DiffusionAssembler_1D.h"

void DiffusionAssembler_1D::computeMassMatrix()
{
    // Assemble mass matrix
	Matrix<double> M(getNumberOfBasisFunctions(), getNumberOfBasisFunctions());
	int N = bspline_x->getKnotvector().distinctKnots.size() - 1; // Number of elements
	for (int ie1 = 0; ie1 < N; ie1++)
	{
		int i_span_1 = bspline_x->getKnotvector().findSpanOfValue(bspline_x->getKnotvector().distinctKnots[ie1]);
		for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
		{
			for (int jl_1 = 0; jl_1 < bspline_x->getDegree() + 1; jl_1++)
			{
				int i1 = i_span_1 - bspline_x->getDegree() + il_1;
				int j1 = i_span_1 - bspline_x->getDegree() + jl_1;

				double v = 0.0;
				std::pair<std::vector<double>, std::vector<double>> gauss = bspline_x->GaussPointsAndWeights(bspline_x->getDegree() + 3, bspline_x->getKnotvector().distinctKnots[ie1], bspline_x->getKnotvector().distinctKnots[ie1 + 1]);
				for (int g1 = 0; g1 < gauss.first.size(); g1++)
				{
					std::pair<std::vector<double>, std::vector<double>> eval = bspline_x->evaluateAtPoint(gauss.first[g1]);
					Matrix<double> J = Jacobian(gauss.first[g1], eval.second);
					double detJ = sqrt(pow(J(0, 0), 2) + pow(J(0, 1), 2));

					double bi_0 = eval.first[il_1];
					double bj_0 = eval.first[jl_1];

					double wvol = gauss.second[g1] * fabs(detJ);

					v += (bi_0 * bj_0) * wvol;
				}

				double temp = M(i1, j1) + v;
				M.setValue(i1, j1, temp);
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
	double start = bspline_x->getKnotvector()(0);
	double end = bspline_x->getKnotvector()(bspline_x->getKnotvector().getSize() - 1);
	double delta = (end - start) / (bspline_x->getNumberOfBasisFunctions() - 1);
	for (int i = 0; i < bspline_x->getNumberOfBasisFunctions() - 1; i++)
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
	Matrix<double> new_sysMatrix(systemMatrix.getNumberOfRows() - boundaryBasisFunctions.size(), systemMatrix.getNumberOfColumns() - boundaryBasisFunctions.size());
	std::vector<double> new_rhs;
	int i = 0;
	for (int ii = 0; ii < stiffnessMatrix.getNumberOfRows(); ii++)
	{
		auto it = std::find_if(boundaryBasisFunctions.begin(), boundaryBasisFunctions.end(), CompareFirst(ii));
		if (it != boundaryBasisFunctions.end())
		{
			continue;
		}
		int j = 0;
		for (int jj = 0; jj < stiffnessMatrix.getNumberOfRows(); jj++)
		{
			it = std::find_if(boundaryBasisFunctions.begin(), boundaryBasisFunctions.end(), CompareFirst(jj));
			if (it != boundaryBasisFunctions.end())
			{
				continue;
			}
			new_sysMatrix.setValue(i, j, systemMatrix(ii, jj));
			j++;
		}
		new_rhs.push_back(rightHandSide[ii]);
		i++;
	}

	systemMatrix = new_sysMatrix;
	rightHandSide = new_rhs;
}

void DiffusionAssembler_1D::applyBoundaryMultipliers()
{
	int new_dim = systemMatrix.getNumberOfRows() + boundaryBasisFunctions.size();
	Matrix<double> new_sysMatrix(new_dim, new_dim, 0.0);
	std::vector<double> new_rhs;
	int cnt2 = 0;
	for (auto it = boundaryBasisFunctions.begin(); it != boundaryBasisFunctions.end(); it++)
	{
		new_sysMatrix.setValue(cnt2, (*it).first + boundaryBasisFunctions.size(), 1.0);
		new_sysMatrix.setValue((*it).first + boundaryBasisFunctions.size(), cnt2, 1.0);
		if (boundaryBasisFunctions[cnt2].second == 1)
		{
			new_rhs.push_back(boundaryConditions->getWestValue());
		}
		else if (boundaryBasisFunctions[cnt2].second == 2)
		{
			new_rhs.push_back(boundaryConditions->getEastValue());
		}
		cnt2++;
	}

	int ii = 0;
	for (int i = cnt2; i < new_sysMatrix.getNumberOfRows(); i++)
	{
		int jj = 0;
		for (int j = cnt2; j < new_sysMatrix.getNumberOfColumns(); j++)
		{
			new_sysMatrix.setValue(i, j, systemMatrix(ii, jj));
			jj++;
		}
		new_rhs.push_back(rightHandSide[ii]);
		ii++;
	}

	systemMatrix = new_sysMatrix;
	rightHandSide = new_rhs;
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