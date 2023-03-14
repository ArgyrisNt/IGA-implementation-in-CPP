#include <iostream>
#include "..\include\Poisson_2D.h"

std::vector<double> Poisson_2D::Ylinspace(int resolution)
{
	std::vector<double> steps;
	double left_limit_y = assembler->getBspline_y().getKnotvector()(0);
	double right_limit_y = assembler->getBspline_y().getKnotvector()(assembler->getBspline_y().getKnotvector().getSize() - 1);
	for (int i = (int)(left_limit_y); i <= resolution; i++)
	{
		double i_step = left_limit_y + (double)(i) * ((right_limit_y - left_limit_y) / ((double)(resolution)));
		steps.push_back(i_step);
	}
	return steps;
}

void Poisson_2D::plotSolution(int resolution)
{
	// Create B-spline surface
	std::vector<double> i_steps = Xlinspace(resolution);
	std::vector<double> j_steps = Ylinspace(resolution);

	std::string filename2("solution.dat");
    std::ofstream my_file2(filename2);
    my_file2 << "variables= " << "\"x\"" << "," << "\"y\"" << "," << "\"sol\"" << "\n";
    my_file2 << "zone t= " << "\"1\"" << ",i=" << resolution+1 << ",j=" << resolution+1 << "\n";

	for (int i = 0; i < i_steps.size(); i++)
	{
		for (int j = 0; j < j_steps.size(); j++)
		{
			if (assembler->trimmingCurve.isCartesianPointInside(i_steps[i], j_steps[j])) continue;

			int span_i = assembler->XspanOfValueInKnotVector(i_steps[i]);
			std::vector<double> bVal_i = assembler->getBspline_x().evaluateAtPoint(i_steps[i]).first;
			int span_j = assembler->YspanOfValueInKnotVector(j_steps[j]);
			std::vector<double> bVal_j = assembler->getBspline_y().evaluateAtPoint(j_steps[j]).first;

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