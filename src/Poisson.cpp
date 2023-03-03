#include <iostream>
#include "..\include\Poisson.h"

template<class T>
void Poisson<T>::constructSol()
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

    std::vector<double> final_sol;
    if (assembler->getBoundaryMode() == "Multipliers")
    {
	    for (int i = 0; i < solution.size(); i++)
	    {
		if (i >= assembler->getBoundaryIds().size())
		        {
			        final_sol.push_back(solution[i]);
		    }
	    }
    }
    else
    {
	    int j = 0;
	    for (int i = 0; i < assembler->getNOF(); i++)
	    {
            auto it = std::find_if(assembler->getBoundaryIds().begin(), assembler->getBoundaryIds().end(), CompareFirst(i));
		    if (it != assembler->getBoundaryIds().end())
		    {	
			    int position = it - assembler->getBoundaryIds().begin();
			    if (assembler->getBoundaryIds()[position].second == 1)
			    {
				    final_sol.push_back(assembler->getBc().getWval());
			    }
			    else if (assembler->getBoundaryIds()[position].second == 2)
			    {
				    final_sol.push_back(assembler->getBc().getEval());
			    }	
                else if (assembler->getBoundaryIds()[position].second == 3)
			    {
				    final_sol.push_back(assembler->getBc().getSval());
			    }
                else if (assembler->getBoundaryIds()[position].second == 4)
			    {
				    final_sol.push_back(assembler->getBc().getNval());
			    }							
		    }
		    else
		    {
			    final_sol.push_back(solution[j]);
			    j++;
		    }
	    }
    }
    solution = final_sol;
}

template<class T>
void Poisson<T>::plotSol(std::string filename1, std::string filename2)
{
    // Create B-spline curve
    double left_limit_x = assembler->getBspline_x().getKnotvector()[0];
    double right_limit_x = assembler->getBspline_x().getKnotvector()[assembler->getBspline_x().getKnotvector().size() - 1];

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
        int span = assembler->getBspline_x().findSpan(i_step);
        std::vector<double> bVal = assembler->getBspline_x().eval(i_step).first;

        double coord_x = 0.0, coord_y = 0.0, coord_z = 0.0;
        for (int kk = 0; kk < bVal.size(); kk++)
        {
            coord_x += bVal[kk] * assembler->getCtrlPts()[span - assembler->getBspline_x().getDegree() + kk][0];
            coord_y += bVal[kk] * assembler->getCtrlPts()[span - assembler->getBspline_x().getDegree() + kk][1];
            coord_z += bVal[kk] * solution[span - assembler->getBspline_x().getDegree() + kk]; // interpolation
        }
        my_file1 << coord_x << " " << coord_y << "\n";
        my_file2 << coord_x << " " << coord_z << "\n";
    }
    my_file1.close();
    my_file2.close();
}