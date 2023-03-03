#include <iostream>
#include "..\include\DiffusionAssembler_1D.h"

void DiffusionAssembler_1D::calcMass()
{
    // Assemble mass matrix
	Matrix<double> M(bspline_x->getNOF(), bspline_x->getNOF());
	int N = bspline_x->knots.size() - 1; // Number of elements
	for (int ie1 = 0; ie1 <= N; ie1++)
	{
		int i_span_1 = bspline_x->findSpan(bspline_x->knots[ie1]);
		for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
		{
			for (int jl_1 = 0; jl_1 < bspline_x->getDegree() + 1; jl_1++)
			{
				int i1 = i_span_1 - bspline_x->getDegree() + il_1;
				int j1 = i_span_1 - bspline_x->getDegree() + jl_1;

				double v = 0.0;
				std::pair<std::vector<double>, std::vector<double>> gauss = bspline_x->calcGaussPts(bspline_x->getDegree() + 3, bspline_x->knots[ie1], bspline_x->knots[ie1 + 1]);
				for (int g1 = 0; g1 < gauss.first.size(); g1++)
				{
					std::pair<std::vector<double>, std::vector<double>> eval = bspline_x->eval(gauss.first[g1]);

					double jacob = calcJacobian(g1, bspline_x->findSpan(gauss.first[g1]), eval.second);

					std::vector<double> ph_bVal;
					for (int kk = 0; kk < eval.first.size(); kk++)
					{
						ph_bVal.push_back((1.0 / jacob) * eval.first[kk]);
					}

					double bi_0 = ph_bVal[il_1];
					double bj_0 = ph_bVal[jl_1];

					double wvol = gauss.second[g1];

					v += jacob * (bi_0 * bj_0) * wvol;
				}

				double temp = M(i1, j1) + v;
				M.setValue(i1, j1, temp);
			}
		}
	}
	mass = M;
}

std::vector<double> DiffusionAssembler_1D::nextStep(std::vector<double> sol)
{
    std::vector<double> b(rhs.size(), 0.0);
	std::vector<double> temp = mass * sol;
	for (int i = 0; i < rhs.size(); i++)
	{
		b[i] = rhs[i] + temp[i];
	}

	return b;
}

std::vector<double> DiffusionAssembler_1D::applyInitCond(double (*func)(double))
{
	std::vector<double> sol;
	std::vector<double> linspaced;
	double start = bspline_x->getKnotvector()[0];
	double end = bspline_x->getKnotvector()[bspline_x->getKnotvector().size() - 1];
	double delta = (end - start) / (bspline_x->getNOF() - 1);
	for (int i = 0; i < bspline_x->getNOF() - 1; i++)
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

void DiffusionAssembler_1D::applyBoundEllimination()
{
    struct CompareFirst
    {
        CompareFirst(int val) : val_(val) {}
        bool operator()(const std::pair<int,char>& elem) const {
            return val_ == elem.first;
        }
    private:
        int val_;
    };

	Matrix<double> new_sysMatrix(sysMat.getRows() - boundary_ids.size(), sysMat.getCols() - boundary_ids.size());
	std::vector<double> new_rhs;
	int i = 0;
	for (int ii = 0; ii < stiff.getRows(); ii++)
	{
		auto it = std::find_if(boundary_ids.begin(), boundary_ids.end(), CompareFirst(ii));
		if (it != boundary_ids.end())
		{
			continue;
		}
		int j = 0;
		for (int jj = 0; jj < stiff.getRows(); jj++)
		{
			it = std::find_if(boundary_ids.begin(), boundary_ids.end(), CompareFirst(jj));
			if (it != boundary_ids.end())
			{
				continue;
			}
			new_sysMatrix.setValue(i, j, sysMat(ii, jj));
			j++;
		}
		new_rhs.push_back(rhs[ii]);
		i++;
	}

	sysMat = new_sysMatrix;
	rhs = new_rhs;
}

void DiffusionAssembler_1D::applyBoundMultipliers()
{
	int new_dim = sysMat.getRows() + boundary_ids.size();
	Matrix<double> new_sysMatrix(new_dim, new_dim, 0.0);
	std::vector<double> new_rhs;
	int cnt2 = 0;
	for (auto it = boundary_ids.begin(); it != boundary_ids.end(); it++)
	{
		new_sysMatrix.setValue(cnt2, (*it).first+boundary_ids.size(), 1.0);
		new_sysMatrix.setValue((*it).first+boundary_ids.size(), cnt2, 1.0);
		if (boundary_ids[cnt2].second == 1)
		{
			new_rhs.push_back(bc->getWval());
		}
		else if (boundary_ids[cnt2].second == 2)
		{
			new_rhs.push_back(bc->getEval());
		}
		cnt2++;
	}

	int ii = 0;
	for (int i = cnt2; i < new_sysMatrix.getRows(); i++)
	{
		int jj = 0;
		for (int j = cnt2; j < new_sysMatrix.getCols(); j++)
		{
			new_sysMatrix.setValue(i, j, sysMat(ii, jj));
			jj++;
		}
		new_rhs.push_back(rhs[ii]);
		ii++;
	}

	sysMat = new_sysMatrix;
	rhs = new_rhs;
}

void DiffusionAssembler_1D::enforceBoundary(std::string& mode)
{
    boundaryMode = mode;
    if (mode == "Ellimination")
	{	
		applyBoundEllimination();
	}
	else if (mode == "Multipliers")
	{
		applyBoundMultipliers();
	}
	else
	{	
		std::cout << "Invalid method for enforcing boundary conditions." << std::endl;
		throw std::invalid_argument("Invalid method");
	}
}