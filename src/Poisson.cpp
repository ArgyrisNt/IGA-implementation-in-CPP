#include <iostream>
#include "..\include\Poisson.h"

template <class T>
void Poisson<T>::expandSolutionOnBoundary()
{
    std::vector<double> newSolution;
    if (assembler->getBoundaryMode() == "Multipliers")
    {
        solution.erase(solution.begin(), solution.begin() + assembler->getBoundaryBasisFunctions().size());
    }
    else
    {
	    int j = 0;
	    for (int i = 0; i < assembler->getNumberOfBasisFunctions(); i++)
	    {
            auto it = std::find_if(assembler->getBoundaryBasisFunctions().begin(), assembler->getBoundaryBasisFunctions().end(), CompareFirst(i));
            if (it != assembler->getBoundaryBasisFunctions().end())
            {
                int position = it - assembler->getBoundaryBasisFunctions().begin();
                if (assembler->getBoundaryBasisFunctions()[position].second == 1)
                {
                    newSolution.push_back(assembler->getBoundaryConditions().getWestValue());
                }
                else if (assembler->getBoundaryBasisFunctions()[position].second == 2)
                {
                    newSolution.push_back(assembler->getBoundaryConditions().getEastValue());
                }
                else if (assembler->getBoundaryBasisFunctions()[position].second == 3)
                {
                    newSolution.push_back(assembler->getBoundaryConditions().getSouthValue());
                }
                else if (assembler->getBoundaryBasisFunctions()[position].second == 4)
                {
                    newSolution.push_back(assembler->getBoundaryConditions().getNorthValue());
                }							
		    }
		    else
		    {
                newSolution.push_back(solution[j]);
                j++;
		    }
	    }
        solution = newSolution;
    }
}

template <class T>
void Poisson<T>::plotSolution(int resolution)
{
    // Create B-spline curve
    double left_limit_x = assembler->getBspline_x().getKnotvector()(0);
    double right_limit_x = assembler->getBspline_x().getKnotvector()(assembler->getBspline_x().getKnotvector().getSize() - 1);

    std::string filename("solution.dat");
    std::ofstream my_file(filename);

    my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
    my_file << "zone t= " << "\"1\"" << ",i=" << resolution + 1 << ",j=" << resolution + 1 << "\n";

    for (int i = (int)(left_limit_x); i <= resolution; i++)
    {
        double i_step = left_limit_x + (double)(i) * ((right_limit_x - left_limit_x) / ((double) (resolution)));
        int span = assembler->XspanOfValueInKnotVector(i_step);
        std::vector<double> bVal = assembler->getBspline_x().evaluateAtPoint(i_step).first;

        double coord_x = 0.0, coord_y = 0.0, coord_z = 0.0;
        for (int kk = 0; kk < bVal.size(); kk++)
        {
            coord_x += bVal[kk] * assembler->getControlPoints()[span - assembler->getBspline_x().getDegree() + kk][0];
            coord_y += bVal[kk] * assembler->getControlPoints()[span - assembler->getBspline_x().getDegree() + kk][1];
            coord_z += bVal[kk] * solution[span - assembler->getBspline_x().getDegree() + kk]; // interpolation
        }
        my_file << coord_x << " " << coord_z << "\n";
    }
    my_file.close();
}