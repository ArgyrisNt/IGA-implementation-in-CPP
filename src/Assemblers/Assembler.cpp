#include <iostream>
#include "..\include\Assembler.h"

Matrix<double> &Assembler::getStiffnessMatrix() 
{ 
	return stiffnessMatrix; 
}

Matrix<double> &Assembler::getSystemMatrix()
{ 
	return systemMatrix;
}

std::vector<double> &Assembler::getRightHandSide()
{
	return rightHandSide;
}

Bspline &Assembler::getBspline_x()
{
	return *bspline_x;
}

double Assembler::getDistinctKnotX(int position)
{
	return bspline_x->getKnotvector().getDistinctKnots()[position];
}

std::vector<double> Assembler::getDistinctKnotsX()
{
	return bspline_x->getKnotvector().getDistinctKnots();
}



int Assembler::XspanOfValueInKnotVector(double value)
{
	return bspline_x->findSpanOfValue(value);
}



void AssemblerBoundary::enforceBoundaryConditions(std::string &mode)
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

void AssemblerBoundary::applyBoundaryEllimination()
{
	int newDimension = systemMatrix.getNumberOfRows() - boundaryBasisFunctions.size();
	Matrix<double> newSystemMatrix(newDimension, newDimension);
	std::vector<double> newRightHandSide;
	int i = 0;
	for (int ii = 0; ii < systemMatrix.getNumberOfRows(); ii++)
	{
		auto it = std::find_if(boundaryBasisFunctions.begin(), boundaryBasisFunctions.end(), CompareFirst(ii));
		if (it != boundaryBasisFunctions.end())
			continue;
		int j = 0;
		for (int jj = 0; jj < systemMatrix.getNumberOfRows(); jj++)
		{
			it = std::find_if(boundaryBasisFunctions.begin(), boundaryBasisFunctions.end(), CompareFirst(jj));
			if (it != boundaryBasisFunctions.end())
				continue;
			newSystemMatrix.setValue(i, j, systemMatrix(ii, jj));
			j++;
		}
		newRightHandSide.push_back(rightHandSide[ii]);
		i++;
	}

	systemMatrix = newSystemMatrix;
	rightHandSide = newRightHandSide;
}

void AssemblerBoundary::applyBoundaryMultipliers()
{
	int new_dim = systemMatrix.getNumberOfRows() + boundaryBasisFunctions.size();
	Matrix<double> newSystemMatrix(new_dim, new_dim);
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

double AssemblerBoundary::addBoundaryValueToRhs(int position)
{
	if (boundaryBasisFunctions[position].second == 1)
	{
		return boundaryConditions->getWestValue();
	}
	else if (boundaryBasisFunctions[position].second == 2)
	{
		return boundaryConditions->getEastValue();
	}
	else if (boundaryBasisFunctions[position].second == 3)
	{
		return boundaryConditions->getNorthValue();
	}
	else if (boundaryBasisFunctions[position].second == 4)
	{
		return boundaryConditions->getSouthValue();
	}
}