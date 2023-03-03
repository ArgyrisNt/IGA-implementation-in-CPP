#include <iostream>
#include "..\include\Assembler.h"


void Assembler::applyBoundEllimination()
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

	Matrix<double> new_stiff(stiff.getRows() - boundary_ids.size(), stiff.getCols() - boundary_ids.size());
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
			new_stiff.setValue(i, j, stiff(ii, jj));
			j++;
		}
		new_rhs.push_back(rhs[ii]);
		i++;
	}

	stiff = new_stiff;
	rhs = new_rhs;
}

void Assembler::applyBoundMultipliers()
{
	int new_dim = stiff.getRows() + boundary_ids.size();
	Matrix<double> new_stiff(new_dim, new_dim);
	std::vector<double> new_rhs;
	int cnt2 = 0;
	for (auto it = boundary_ids.begin(); it != boundary_ids.end(); it++)
	{
		new_stiff.setValue(cnt2, (*it).first+boundary_ids.size(), 1.0);
		new_stiff.setValue((*it).first+boundary_ids.size(), cnt2, 1.0);
		if (boundary_ids[cnt2].second == 1)
		{
			new_rhs.push_back(bc->getWval());
		}
		else if (boundary_ids[cnt2].second == 2)
		{
			new_rhs.push_back(bc->getEval());
		}
		else if (boundary_ids[cnt2].second == 3)
		{
			new_rhs.push_back(bc->getNval());
		}
		else if (boundary_ids[cnt2].second == 4)
		{
			new_rhs.push_back(bc->getSval());
		}
		cnt2++;
	}

	int ii = 0;
	for (int i = cnt2; i < new_stiff.getRows(); i++)
	{
		int jj = 0;
		for (int j = cnt2; j < new_stiff.getCols(); j++)
		{
			new_stiff.setValue(i, j, stiff(ii, jj));
			jj++;
		}
		new_rhs.push_back(rhs[ii]);
		ii++;
	}

	stiff = new_stiff;
	rhs = new_rhs;
}

void Assembler::enforceBoundary(std::string& mode)
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