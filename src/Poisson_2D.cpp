#include <iostream>
#include "..\include\Poisson_2D.h"

void Poisson_2D::plotSolution(int resolution)
{
	// Create B-spline surface
	double left_limit_x = assembler->getBspline_x().getKnotvector()(0);
	double right_limit_x = assembler->getBspline_x().getKnotvector()(assembler->getBspline_x().getKnotvector().getSize() - 1);
	double left_limit_y = assembler->getBspline_y().getKnotvector()(0);
	double right_limit_y = assembler->getBspline_y().getKnotvector()(assembler->getBspline_y().getKnotvector().getSize() - 1);

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

			bool isInsideTrimmingCurve = (std::pow(i_step - assembler->trimmingCurve.center.x, 2) + std::pow(j_step - assembler->trimmingCurve.center.y, 2)) < std::pow(assembler->trimmingCurve.radius, 2);
			if (isInsideTrimmingCurve) continue;

			int span_i = assembler->getBspline_x().getKnotvector().findSpanOfValue(i_step);
			std::vector<double> bVal_i = assembler->getBspline_x().evaluateAtPoint(i_step).first;
			int span_j = assembler->getBspline_y().getKnotvector().findSpanOfValue(j_step);
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
			my_file2 << coord_x << " " << coord_y << " " << coord_z << "\n";
		}
	}
	my_file2.close();

	std::string filename3("parameter_space.dat");
	assembler->writeParameterSpaceToFile(filename3);

	if (assembler->trimmed_triangles.size())
	{
		std::string filename4("trimmed_triangles.obj");
		writeTrimmedTrianglesToFile(filename4);
	}
}

void Poisson_2D::writeTrimmedTrianglesToFile(std::string filename)
{
	std::ofstream my_file(filename);
	my_file << "# " << assembler->trimmed_triangles.size() << "\n";
	int cnt = 0;
	for (auto triangle : assembler->trimmed_triangles)
	{
		cnt += 2;
		my_file << "v " << triangle.vertex1.x << " " << triangle.vertex1.y << "\n";
		my_file << "v " << triangle.vertex2.x << " " << triangle.vertex2.y << "\n";
		my_file << "l " << cnt - 1 << " " << cnt << "\n";
		cnt += 2;
		my_file << "v " << triangle.vertex2.x << " " << triangle.vertex2.y << "\n";
		my_file << "v " << triangle.vertex3.x << " " << triangle.vertex3.y << "\n";
		my_file << "l " << cnt - 1 << " " << cnt << "\n";
		cnt += 2;
		my_file << "v " << triangle.vertex3.x << " " << triangle.vertex3.y << "\n";
		my_file << "v " << triangle.vertex1.x << " " << triangle.vertex1.y << "\n";
		my_file << "l " << cnt - 1 << " " << cnt << "\n";
	}
	my_file.close();
}