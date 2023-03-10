#include <iostream>
#include "..\include\Assembler.h"

void Assembler::XcomputeDistinctKnots()
{
	XdistinctKnots = {};
	double currentValue, previousValue = -100.0;
	KnotVector<double> knotVector = bspline_x->getKnotvector();
	for (int i = bspline_x->getDegree(); i < knotVector.getSize() - bspline_x->getDegree(); i++)
	{
		currentValue = knotVector(i);
		if (currentValue != previousValue)
			XdistinctKnots.push_back(knotVector(i));
		previousValue = knotVector(i);
	}
}

int Assembler::XspanOfValueInKnotVector(double value)
{
	return bspline_x->getKnotvector().findSpanOfValue(value);
}

void Assembler::applyBoundaryEllimination()
{
	int newDimension = stiffnessMatrix.getNumberOfRows() - boundaryBasisFunctions.size();
	Matrix<double> newStiffnessMatrix(newDimension, newDimension);
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

double Assembler::addBoundaryValueToRhs(int position)
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
		newRightHandSide.push_back(addBoundaryValueToRhs(cnt2));		
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

void Assembler::mapValuesToDomain(std::vector<double> &GaussPoints, const double left, const double right)
{
	for (int i = 0; i < GaussPoints.size(); i++)
	{
		GaussPoints[i] = ((left * (1 - GaussPoints[i]) + right * (1 + GaussPoints[i])) / 2);
	}
}

std::vector<std::pair<double, double>> Assembler::GaussPointsAndWeights(int numberOfPoints, const double left, const double right)
{
	assert(numberOfPoints > 0);
	if (numberOfPoints > 5) numberOfPoints = 5;
	std::vector<double> GaussPoints, GaussWeights;

	// Compute Gauss points in interval [0,1]
	switch (numberOfPoints)
	{
	case 1:
		GaussPoints = {0.0};
		GaussWeights = {2.0};
		break;
	case 2:
		GaussPoints = {-0.57735, 0.57735};
		GaussWeights = {1.0, 1.0};
		break;
	case 3:
		GaussPoints = {0.0, -0.774597, 0.774597};
		GaussWeights = {0.888889, 0.555556, 0.555556};
		break;
	case 4:
		GaussPoints = {-0.861136, -0.339981, 0.339981, 0.861136};
		GaussWeights = {0.347855, 0.652145, 0.652145, 0.347855};
		break;
	case 5:
		GaussPoints = {-0.90618, -0.538469, 0.0, 0.538469, 0.90618};
		GaussWeights = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927};
		break;
	default:
		std::cout << "Invalid dimension. Valid dimensions are 1,2,3,4,5." << std::endl;
		throw std::invalid_argument("Invalid dimension");
		break;
	}
	mapValuesToDomain(GaussPoints, left, right);

	std::vector<std::pair<double, double>> GaussPointsAndWeights;
	for (int i = 0; i < GaussPoints.size(); i++)
	{
		GaussPointsAndWeights.push_back(std::make_pair(GaussPoints[i], GaussWeights[i]));
	}

	return GaussPointsAndWeights;
}