#include <iostream>
#include "..\include\Poisson_2D.h"

void Poisson_2D::plotSolution(int resolution)
{
	// Create B-spline surface
	double left_limit_x = assembler->getBspline_x().getKnotvector()[0];
	double right_limit_x = assembler->getBspline_x().getKnotvector()[assembler->getBspline_x().getKnotvector().size() - 1];
	double left_limit_y = assembler->getBspline_y().getKnotvector()[0];
	double right_limit_y = assembler->getBspline_y().getKnotvector()[assembler->getBspline_y().getKnotvector().size() - 1];

	std::string filename1("surface.dat");
    std::ofstream my_file1(filename1);
    my_file1 << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
    my_file1 << "zone t= " << "\"1\"" << ",i=" << resolution+1 << ",j=" << resolution+1 << "\n";

	std::string filename2("solution.dat");
    std::ofstream my_file2(filename2);
    my_file2 << "variables= " << "\"x\"" << "," << "\"y\"" << "," << "\"sol\"" << "\n";
    my_file2 << "zone t= " << "\"1\"" << ",i=" << resolution+1 << ",j=" << resolution+1 << "\n";

	for (int i = (int)(left_limit_x); i <= resolution; i++)
	{
		for (int j = (int)(left_limit_x); j <= resolution; j++)
		{
			double i_step = left_limit_x + (double)(i) * ((right_limit_x - left_limit_x) / ((double) (resolution)));
			double j_step = left_limit_y + (double)(j) * ((right_limit_y - left_limit_y) / ((double) (resolution)));

			if ((std::pow(i_step - assembler->trimming[0], 2) + std::pow(j_step - assembler->trimming[1], 2)) < std::pow(assembler->trimming[2], 2)) continue;

			int span_i = assembler->getBspline_x().findSpanInVector(i_step);
			std::vector<double> bVal_i = assembler->getBspline_x().evaluateAtPoint(i_step).first;
			int span_j = assembler->getBspline_y().findSpanInVector(j_step);
			std::vector<double> bVal_j = assembler->getBspline_y().evaluateAtPoint(j_step).first;

			double coord_x = 0.0, coord_y = 0.0, coord_z = 0.0;
			for (int kkx = 0; kkx < bVal_i.size(); kkx++)
			{
				for (int kky = 0; kky < bVal_j.size(); kky++)
				{
					int i1 = span_i - assembler->getBspline_x().getDegree() + kkx;
					int i2 = span_j - assembler->getBspline_y().getDegree() + kky;
					int my = i1 * assembler->getBspline_y().getNumberOfBasisFunctions() + i2;
					coord_x += bVal_i[kkx] * bVal_j[kky] * assembler->getControlPoints()[my][0];
					coord_y += bVal_i[kkx] * bVal_j[kky] * assembler->getControlPoints()[my][1];
					coord_z += bVal_i[kkx] * bVal_j[kky] * solution[my];
				}
			}
			my_file1 << coord_x << " " << coord_y << "\n";
			my_file2 << coord_x << " " << coord_y << " " << coord_z << "\n";
		}
	}
	my_file1.close();
	my_file2.close();

	std::string filename3("parameter_space.dat");
	std::ofstream my_file3(filename3);
	my_file3 << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
	my_file3 << "zone t= " << "\"1\"" << ",i=" << assembler->getBspline_x().distinctKnots.size() << ",j=" << assembler->getBspline_y().distinctKnots.size() << "\n";
	for (int j = 0; j < assembler->getBspline_y().distinctKnots.size(); j++)
	{
		for (int i = 0; i < assembler->getBspline_x().distinctKnots.size(); i++)
		{
			my_file3 << assembler->getBspline_x().distinctKnots[i] << " " << assembler->getBspline_y().distinctKnots[j] << "\n";
		}
	}
	my_file3.close();

	std::string filename4("trimmed_triangles.obj");
	std::ofstream my_file4(filename4);

	my_file4 << "# " << assembler->trimmed_triangles.size() << "\n";
	int cnt = 0;
	for (auto triangle : assembler->trimmed_triangles)
	{
		for (int i = 0; i < 3; i ++)
		{
			cnt += 2;
			my_file4 << "v " << triangle[i].first << " " << triangle[i].second << "\n";
			my_file4 << "v " << triangle[(i + 1) % 3].first << " " << triangle[(i + 1) % 3].second << "\n";
			my_file4 << "l " << cnt - 1 << " " << cnt << "\n";
		}
	}
	my_file4.close();
}