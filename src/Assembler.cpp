#include <iostream>
#include "..\include\Assembler.h"

template <class T>
Matrix<double> &Assembler<T>::getStiffnessMatrix()
{ 
	return stiffnessMatrix; 
}

template <class T>
Matrix<double> &Assembler<T>::getSystemMatrix()
{ 
	return systemMatrix;
}

template <class T>
std::vector<double> &Assembler<T>::getRightHandSide()
{
	return rightHandSide;
}

template <class T>
Bspline &Assembler<T>::getBsplineEntity()
{
	return bsplineEntity;
}

template <class T>
double Assembler<T>::getDistinctKnotX(int position)
{
	return bsplineEntity.getKnotvector().getDistinctKnots()[position];
}

template <class T>
std::vector<double> Assembler<T>::getDistinctKnotsX()
{
	return bsplineEntity.getKnotvector().getDistinctKnots();
}



template <class T>
int Assembler<T>::XspanOfValueInKnotVector(double value)
{
	return bsplineEntity.findSpanOfValue(value);
}



template <class T>
void Assembler<T>::enforceBoundaryConditions(std::string &mode)
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

template <class T>
void Assembler<T>::applyBoundaryEllimination()
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

template <class T>
void Assembler<T>::applyBoundaryMultipliers()
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

template <class T>
double Assembler<T>::addBoundaryValueToRhs(int position)
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