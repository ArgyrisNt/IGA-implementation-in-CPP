#include <iostream>
#include "..\include\Assembler.h"


void Assembler::applyBoundaryEllimination()
{
	Matrix<double> newStiffnessMatrix(stiffnessMatrix.getNumberOfRows() - boundaryBasisFunctions.size(), stiffnessMatrix.getNumberOfColumns() - boundaryBasisFunctions.size());
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
			newStiffnessMatrix.setValue(i, j, stiffnessMatrix(ii, jj));
			j++;
		}
		newRightHandSide.push_back(rightHandSide[ii]);
		i++;
	}

	stiffnessMatrix = newStiffnessMatrix;
	rightHandSide = newRightHandSide;
}

void Assembler::applyBoundaryMultipliers()
{
	int new_dim = stiffnessMatrix.getNumberOfRows() + boundaryBasisFunctions.size();
	Matrix<double> newStiffnessMatrix(new_dim, new_dim);
	std::vector<double> newRightHandSide;
	int cnt2 = 0;
	for (auto it = boundaryBasisFunctions.begin(); it != boundaryBasisFunctions.end(); it++)
	{
		newStiffnessMatrix.setValue(cnt2, (*it).first + boundaryBasisFunctions.size(), 1.0);
		newStiffnessMatrix.setValue((*it).first + boundaryBasisFunctions.size(), cnt2, 1.0);
		if (boundaryBasisFunctions[cnt2].second == 1)
		{
			newRightHandSide.push_back(boundaryConditions->getWestValue());
		}
		else if (boundaryBasisFunctions[cnt2].second == 2)
		{
			newRightHandSide.push_back(boundaryConditions->getEastValue());
		}
		else if (boundaryBasisFunctions[cnt2].second == 3)
		{
			newRightHandSide.push_back(boundaryConditions->getNorthValue());
		}
		else if (boundaryBasisFunctions[cnt2].second == 4)
		{
			newRightHandSide.push_back(boundaryConditions->getSouthValue());
		}
		cnt2++;
	}

	int ii = 0;
	for (int i = cnt2; i < newStiffnessMatrix.getNumberOfRows(); i++)
	{
		int jj = 0;
		for (int j = cnt2; j < newStiffnessMatrix.getNumberOfColumns(); j++)
		{
			newStiffnessMatrix.setValue(i, j, stiffnessMatrix(ii, jj));
			jj++;
		}
		newRightHandSide.push_back(rightHandSide[ii]);
		ii++;
	}

	stiffnessMatrix = newStiffnessMatrix;
	rightHandSide = newRightHandSide;
}

void Assembler::enforceBoundaryConditions(std::string& mode)
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