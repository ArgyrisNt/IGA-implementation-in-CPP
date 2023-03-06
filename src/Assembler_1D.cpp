#include <iostream>
#include <cmath>
#include "..\include\Assembler_1D.h"

void Assembler_1D::computeBoundary()
{
	for (int i = 0; i < getNumberOfBasisFunctions(); i++)
	{
		if ((i == 0) && (boundaryConditions->getWestType() == "Dirichlet"))
		{
			boundaryBasisFunctions.push_back(std::make_pair(i, 1));
		}
		else if ((i == getNumberOfBasisFunctions() - 1) && (boundaryConditions->getEastType() == "Dirichlet"))
		{
			boundaryBasisFunctions.push_back(std::make_pair(i, 2));
		}
	}
}

double Assembler_1D::Jacobian(double g1, int span_g1, std::vector<double> &Nx)
{	
	double temp1 = 0.0, temp2 = 0.0;
	for (int kk = 0; kk < Nx.size(); kk++)
	{
		temp1 += Nx[kk] * getControlPoints()[span_g1 - bspline_x->getDegree() + kk][0];
		temp2 += Nx[kk] * getControlPoints()[span_g1 - bspline_x->getDegree() + kk][1];
	}
	double jacob1 = sqrt(pow(temp1,2) + pow(temp2,2));
	
	return jacob1;
}

void Assembler_1D::computeStiffnessMatrix()
{
	// Assemble stiffness matrix
	Matrix<double> A(getNumberOfBasisFunctions(), getNumberOfBasisFunctions());
	int N = bspline_x->distinctKnots.size() - 1; // Number of elements
	for (int ie1 = 0; ie1 < N; ie1++)
	{
		int i_span_1 = bspline_x->findSpanInVector(bspline_x->distinctKnots[ie1]);
		for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
		{
			for (int jl_1 = 0; jl_1 < bspline_x->getDegree() + 1; jl_1++)
			{
				int i1 = i_span_1 - bspline_x->getDegree() + il_1;
				int j1 = i_span_1 - bspline_x->getDegree() + jl_1;

				double v = 0.0;
				std::pair<std::vector<double>, std::vector<double>> gauss = bspline_x->GaussPointsAndWeights(bspline_x->getDegree() + 3, bspline_x->distinctKnots[ie1], bspline_x->distinctKnots[ie1 + 1]);
				for (int g1 = 0; g1 < gauss.first.size(); g1++)
				{
					std::vector<double> gradVal = bspline_x->evaluateAtPoint(gauss.first[g1]).second;
					double jacob = Jacobian(g1, bspline_x->findSpanInVector(gauss.first[g1]), gradVal);

					std::vector<double> ph_grad_val;
					for (int kk = 0; kk < gradVal.size(); kk++)
					{
						ph_grad_val.push_back((1.0 / jacob) * gradVal[kk]);
					}

					double bi_x = ph_grad_val[il_1];
					double bj_x = ph_grad_val[jl_1];

					double wvol = gauss.second[g1];

					v += jacob * (bi_x * bj_x) * wvol;
				}

				double temp = A(i1, j1) + v;
				A.setValue(i1, j1, temp);
			}
		}
	}
	stiffnessMatrix = A;
}

void Assembler_1D::computeRightHandSide()
{
	// Assemble rhs vector
	std::vector<double> b(getNumberOfBasisFunctions(), 0.0);
	int N = bspline_x->distinctKnots.size() - 1;
	for (int ie1 = 0; ie1 < N; ie1++)
	{
		int i_span_1 = bspline_x->findSpanInVector(bspline_x->distinctKnots[ie1]);
		for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
		{
			int i1 = i_span_1 - bspline_x->getDegree() + il_1;

			double v = 0.0;
			std::pair<std::vector<double>, std::vector<double>> gauss = bspline_x->GaussPointsAndWeights(bspline_x->getDegree() + 3, bspline_x->distinctKnots[ie1], bspline_x->distinctKnots[ie1 + 1]);
			for (int g1 = 0; g1 < gauss.first.size(); g1++)
			{
				std::pair<std::vector<double>, std::vector<double>> eval = bspline_x->evaluateAtPoint(gauss.first[g1]);

				double jacob = Jacobian(g1, bspline_x->findSpanInVector(gauss.first[g1]), eval.second);
					
				std::vector<double> ph_val;
				for (int kk = 0; kk < eval.first.size(); kk++)
				{
					ph_val.push_back((1.0 / jacob) * eval.first[kk]);
				}

				double bi_0 = ph_val[il_1]; //ph_bVal[il_1];
				double wvol = gauss.second[g1];

				v += jacob * bi_0 * sourceFunction * wvol;
			}

			b[i1] += v;
		}
	}
	rightHandSide = b;
}
