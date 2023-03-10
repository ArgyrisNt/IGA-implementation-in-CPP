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
void Poisson<T>::plotSolution(std::string filename1, std::string filename2)
{
    // Create B-spline curve
    double left_limit_x = assembler->getBspline_x().getKnotvector()(0);
    double right_limit_x = assembler->getBspline_x().getKnotvector()(assembler->getBspline_x().getKnotvector().getSize() - 1);

    //std::string filename1("curve.dat");
    std::ofstream my_file1(filename1);

    //std::string filename2("solution.dat");
    std::ofstream my_file2(filename2);

    my_file1 << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
    my_file1 << "zone t= " << "\"1\"" << ",i=" << 101 << ",j=" << 101 << "\n";
    my_file2 << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
    my_file2 << "zone t= " << "\"1\"" << ",i=" << 101 << ",j=" << 101 << "\n";

    for (int i = (int) (left_limit_x); i <= 100; i++)
    {
        double i_step = left_limit_x + (double)(i) * ((right_limit_x - left_limit_x) / 100.0);
        int span = assembler->getBspline_x().getKnotvector().findSpanOfValue(i_step);
        std::vector<double> bVal = assembler->getBspline_x().evaluateAtPoint(i_step).first;

        double coord_x = 0.0, coord_y = 0.0, coord_z = 0.0;
        for (int kk = 0; kk < bVal.size(); kk++)
        {
            coord_x += bVal[kk] * assembler->getControlPoints()[span - assembler->getBspline_x().getDegree() + kk][0];
            coord_y += bVal[kk] * assembler->getControlPoints()[span - assembler->getBspline_x().getDegree() + kk][1];
            coord_z += bVal[kk] * solution[span - assembler->getBspline_x().getDegree() + kk]; // interpolation
        }
        my_file1 << coord_x << " " << coord_y << "\n";
        my_file2 << coord_x << " " << coord_z << "\n";
    }
    my_file1.close();
    my_file2.close();
}