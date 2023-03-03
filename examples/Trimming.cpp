#include <iostream>
#include "..\IGA.h"

bool trimmed(double x, double y) { return (std::pow(x - 1.0, 2) + std::pow(y, 2)) < std::pow(0.5, 2); }

int main()
{
    int p = 2;
    std::vector<double> U{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    std::vector<double> W{1.0, 1.0, 1.0};
    Bspline bspline_x(p, U, W);
    Bspline bspline_y(bspline_x);
    std::vector<std::vector<double>> ctrlPts{ {0.0, 0.0}, {0.0, 2.0}, {0.0, 4.0},
                                              {2.0, 0.0}, {2.0, 2.0}, {2.0, 4.0},
                                              {4.0, 0.0}, {4.0, 2.0}, {4.0, 4.0}};
    
    std::vector<std::vector<double>> C;
    double left_limit_x = bspline_x.getKnotvector()[0];
    double right_limit_x = bspline_x.getKnotvector()[bspline_x.getKnotvector().size() - 1];
    double left_limit_y = bspline_y.getKnotvector()[0];
    double right_limit_y = bspline_y.getKnotvector()[bspline_y.getKnotvector().size() - 1];
    std::string filename("trimming.dat");
    std::ofstream my_file(filename);
    my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
    my_file << "zone t= " << "\"1\"" << ",i=" << 101 << ",j=" << 101 << "\n";
    for (int j = (int) (left_limit_y); j <= 100; j++)
    {
        for (int i = (int) (left_limit_x); i <= 100; i++)
        {
            double i_step = left_limit_x + (double) (i) * ((right_limit_x - left_limit_x) / 100.0);
            double j_step = left_limit_y + (double) (j) * ((right_limit_y - left_limit_y) / 100.0);
            if (trimmed(i_step, j_step)) continue;
            int span_i = bspline_x.findSpan(i_step);
            int span_j = bspline_y.findSpan(j_step);
            std::vector<double> bVal_i = bspline_x.eval(i_step).first;
            std::vector<double> bVal_j = bspline_y.eval(j_step).first;

            double coord_x = 0.0, coord_y = 0.0;
            for (int ii = 0; ii < bVal_i.size(); ii++)
			{
				for (int jj = 0; jj < bVal_j.size(); jj++)
				{
                    int i1 = span_i - bspline_x.getDegree() + ii;
                    int j1 = span_j - bspline_y.getDegree() + jj;
                    int index = i1 * bspline_y.getNOF() + j1;
                    coord_x += bVal_i[ii] * bVal_j[jj] * ctrlPts[index][0];
                    coord_y += bVal_i[ii] * bVal_j[jj] * ctrlPts[index][1];
                }
            }
            C.push_back({coord_x, coord_y});
            my_file << coord_x << " " << coord_y << "\n";
        }
    }
    my_file.close();

    return 0;
}