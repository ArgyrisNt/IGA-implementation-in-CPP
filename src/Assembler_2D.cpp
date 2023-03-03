#include <iostream>
#include "..\include\Assembler_2D.h"

struct CompareFirst
{
	CompareFirst(int val) : val_(val) {}
	bool operator()(const std::pair<int, char> &elem) const
	{
		return val_ == elem.first;
	}

private:
	int val_;
};

void Assembler_2D::calcBound()
{
	for (int i = 0; i < getNOF(); i++)
	{
		if ((i >= 0 && i <= bspline_y->getNOF() - 1) && (bc->getWest() == "Dirichlet"))
		{
			auto it = std::find_if(boundary_ids.begin(), boundary_ids.end(), CompareFirst(i));
			if (it == boundary_ids.end())
			{
				boundary_ids.push_back(std::make_pair(i, 1));
			}
		}
		else if ((i % bspline_y->getNOF() == 0) && (bc->getSouth() == "Dirichlet"))
		{
			auto it = std::find_if(boundary_ids.begin(), boundary_ids.end(), CompareFirst(i));
			if (it == boundary_ids.end())
			{
				boundary_ids.push_back(std::make_pair(i,4));
			}
		}
		else if ((i % bspline_y->getNOF() == bspline_y->getNOF() - 1) && (bc->getNorth() == "Dirichlet"))
		{
			auto it = std::find_if(boundary_ids.begin(), boundary_ids.end(), CompareFirst(i));
			if (it == boundary_ids.end())
			{
				boundary_ids.push_back(std::make_pair(i,3));
			}
		}
		else if ((i >= getNOF() - bspline_y->getNOF()) && (bc->getEast() == "Dirichlet"))
		{
			auto it = std::find_if(boundary_ids.begin(), boundary_ids.end(), CompareFirst(i));
			if (it == boundary_ids.end())
			{
				boundary_ids.push_back(std::make_pair(i,2));
			}
		}
	}
}

std::vector<double> Assembler_2D::createTensorProduct(std::vector<double>& vec1, std::vector<double>& vec2)
{
	std::vector<double> tensor_product;
	for (int u = 0; u < vec1.size(); u++)
	{
		for (int v = 0; v < vec2.size(); v++)
		{
			tensor_product.push_back({ vec1[u] * vec2[v] });
		}
	}

	return tensor_product;
}

Matrix<double> Assembler_2D::calcJacobian(double g1, double g2, int span_g1, int span_g2, std::vector<double>& shp_fnc_dx, std::vector<double>& shp_fnc_dy)
{	
	std::vector<int> index, index_x, index_y;
	for (int kk = 0; kk < bspline_x->getDegree() + 1; kk++)
	{
		index_x.push_back(span_g1 - bspline_x->getDegree() + kk);
	}
	for (int kk = 0; kk < bspline_y->getDegree() + 1; kk++)
	{
		index_y.push_back(span_g2 - bspline_y->getDegree() + kk);
	}
	for (auto indx: index_x)
	{
		for (auto indy: index_y)
		{
			index.push_back(indx * (bspline_y->getKnotvector().size() - bspline_y->getDegree() - 1) + indy);
		}
	}
	
	std::vector<double>	temp_x(2, 0.0);
	for (int kk = 0; kk < shp_fnc_dx.size(); kk++)
	{
		temp_x[0] += shp_fnc_dx[kk] * ctrlPts[index[kk]][0];
		temp_x[1] += shp_fnc_dx[kk] * ctrlPts[index[kk]][1];
	}

	std::vector<double>	temp_y(2, 0.0);
	for (int kk = 0; kk < shp_fnc_dy.size(); kk++)
	{
		temp_y[0] += shp_fnc_dy[kk] * ctrlPts[index[kk]][0];
		temp_y[1] += shp_fnc_dy[kk] * ctrlPts[index[kk]][1];
	}

	Matrix<double> Jacobian(2, 2);
	Jacobian.setValue(0, 0, temp_x[0]);
	Jacobian.setValue(0, 1, temp_y[0]);
	Jacobian.setValue(1, 0, temp_x[1]);
	Jacobian.setValue(1, 1, temp_y[1]);
	
	return Jacobian;
}

std::pair<std::vector<double>,std::vector<double>> Assembler_2D::Map2Physical(Matrix<double>& Jacobian, std::vector<double>& shp_fnc_dx, std::vector<double>& shp_fnc_dy)
{
	Matrix<double> inv_Jacobian = Jacobian.inverse();
	std::vector<double> ph_grad_val_1(shp_fnc_dx.size(), 0.0);
	std::vector<double> ph_grad_val_2(shp_fnc_dy.size(), 0.0);
	for (int kk = 0; kk < shp_fnc_dx.size(); kk++)
	{
		std::vector<double> res(2, 0.0);
		std::vector<double> vec{ shp_fnc_dx[kk], shp_fnc_dy[kk] };
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				res[i] +=  vec[j] * inv_Jacobian(i, j);
			}
		}
		ph_grad_val_1[kk] = res[0];
		ph_grad_val_2[kk] = res[1];
	}

	return std::make_pair(ph_grad_val_1, ph_grad_val_2);
}

void Assembler_2D::calcTrimmed()
{
	int Nx = bspline_x->knots.size() - 1; // Number of elements on x-direction
	int Ny = bspline_y->knots.size() - 1; // NUmber of elements on y-direction
	std::vector<std::vector<std::pair<double, double>>> temp_vec(Nx * Ny);
	trimmed_elements_info = temp_vec;
	elements_vertices = temp_vec;
	int elem = -1;
	for (int ie1 = 0; ie1 < Nx; ie1++)
	{
		for (int ie2 = 0; ie2 < Ny; ie2++)
		{
			elem++;
			if (trimming[2] == 0.0)
			{
				trimmed_elements.push_back(std::make_pair(false, 0));
				continue;
			}
			std::pair<double, double> center = std::make_pair((bspline_x->knots[ie1 + 1] + bspline_x->knots[ie1]) / 2.0, (bspline_y->knots[ie2 + 1] + bspline_y->knots[ie2]) / 2.0);
			double di = projection_on_trimming(center);
			double r_in = sqrt(std::pow(bspline_x->knots[ie1 + 1] - center.first, 2));
			double r_out = sqrt(std::pow(bspline_x->knots[ie1] - center.first, 2) + std::pow(bspline_y->knots[ie2] - center.second, 2));

			int count = 0;
			if (checkIfOutside(std::make_pair(bspline_x->knots[ie1], bspline_y->knots[ie2])))
			{
				count++;
			}
			else trimmed_elements_info[elem].push_back(std::make_pair(bspline_x->knots[ie1], bspline_y->knots[ie2]));
			elements_vertices[elem].push_back(std::make_pair(bspline_x->knots[ie1], bspline_y->knots[ie2]));
			if (checkIfOutside(std::make_pair(bspline_x->knots[ie1], bspline_y->knots[ie2 + 1])))
			{
				count++;
			}
			else trimmed_elements_info[elem].push_back(std::make_pair(bspline_x->knots[ie1], bspline_y->knots[ie2 + 1]));
			elements_vertices[elem].push_back(std::make_pair(bspline_x->knots[ie1], bspline_y->knots[ie2 + 1]));
			if (checkIfOutside(std::make_pair(bspline_x->knots[ie1 + 1], bspline_y->knots[ie2]))) 
			{
				count++;
			}
			else trimmed_elements_info[elem].push_back(std::make_pair(bspline_x->knots[ie1 + 1], bspline_y->knots[ie2]));
			elements_vertices[elem].push_back(std::make_pair(bspline_x->knots[ie1 + 1], bspline_y->knots[ie2]));
			if (checkIfOutside(std::make_pair(bspline_x->knots[ie1 + 1], bspline_y->knots[ie2 + 1])))
			{
				count++;
			}
			else trimmed_elements_info[elem].push_back(std::make_pair(bspline_x->knots[ie1 + 1], bspline_y->knots[ie2 + 1]));
			elements_vertices[elem].push_back(std::make_pair(bspline_x->knots[ie1 + 1], bspline_y->knots[ie2 + 1]));

			if (di < r_in)
			{			
				trimmed_elements.push_back(std::make_pair(true, count));
			}
			else if (di > r_out) 
			{
				trimmed_elements.push_back(std::make_pair(false, count));
			}
			else if (di > r_in && di < r_out)
			{
				if (count != 0) trimmed_elements.push_back(std::make_pair(true, count));
				else trimmed_elements.push_back(std::make_pair(false, count));
			}
		}
	}
}

void Assembler_2D::calcStiff()
{
	// Assemble stiffness matrix
	int Nx = bspline_x->knots.size() - 1; // Number of elements on x-direction
	int Ny = bspline_y->knots.size() - 1; // NUmber of elements on y-direction
	Matrix<double> A(getNOF(), getNOF());
	int cnt = -1;
	for (int ie1 = 0; ie1 < Nx; ie1++)
	{
		int i_span_1 = bspline_x->findSpan(bspline_x->knots[ie1]);
		for (int ie2 = 0; ie2 < Ny; ie2++)
		{
			cnt++;
			if (trimmed_elements[cnt].first)
			{
				std::vector<std::vector<std::pair<double, double>>> triangles;
				if (trimmed_elements[cnt].second == 1)
				{
					triangles = construct_3_triangles(cnt);
				}
				else if (trimmed_elements[cnt].second == 2)
				{
					triangles = construct_2_triangles(cnt);
				}
				else if (trimmed_elements[cnt].second == 3)
				{
					triangles = construct_1_triangle(cnt);
				}
				for (auto triangle : triangles)
				{
					trimmed_triangles.push_back(triangle);
/*----SOS---->*/	double a_x = std::min({triangle[0].first, triangle[1].first, triangle[2].first});
					double b_x = std::max({triangle[0].first, triangle[1].first, triangle[2].first});
					double a_y = std::min({triangle[0].second, triangle[1].second, triangle[2].second});
					double b_y = std::max({triangle[0].second, triangle[1].second, triangle[2].second});
					int i_span_2 = bspline_y->findSpan(bspline_y->knots[ie2]);
					int mycnt_1 = 0;
					for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
					{
						for (int il_2 = 0; il_2 < bspline_y->getDegree() + 1; il_2++)
						{
							int mycnt_2 = 0;
							for (int jl_1 = 0; jl_1 < bspline_x->getDegree() + 1; jl_1++)
							{
								for (int jl_2 = 0; jl_2 < bspline_y->getDegree() + 1; jl_2++)
								{
									int i1 = i_span_1 - bspline_x->getDegree() + il_1;
									int j1 = i_span_1 - bspline_x->getDegree() + jl_1;

									int i2 = i_span_2 - bspline_y->getDegree() + il_2;
									int j2 = i_span_2 - bspline_y->getDegree() + jl_2;

									double v = 0.0;
									std::pair<std::vector<double>, std::vector<double>> gauss_x = calcGaussPtsTria(a_x, b_x);
									std::pair<std::vector<double>, std::vector<double>> gauss_y = calcGaussPtsTria(a_y, b_y);
									
									for (int g1 = 0; g1 < gauss_x.first.size(); g1++)
									{
										for (int g2 = 0; g2 < gauss_y.first.size(); g2++)
										{
											if (g1 == (gauss_x.first.size()-1) && g2 == (gauss_y.first.size()-1)) break;
											std::pair<std::vector<double>, std::vector<double>> eval_1 = bspline_x->eval(gauss_x.first[g1]);
											std::pair<std::vector<double>, std::vector<double>> eval_2 = bspline_y->eval(gauss_y.first[g2]);

											std::vector<double> shp_fnc_dx = createTensorProduct(eval_1.second, eval_2.first);
											std::vector<double> shp_fnc_dy = createTensorProduct(eval_1.first, eval_2.second);

											Matrix<double> Jacobian = calcJacobian(g1, g2, bspline_x->findSpan(gauss_x.first[g1]), bspline_y->findSpan(gauss_y.first[g2]), shp_fnc_dx, shp_fnc_dy);
											std::pair<std::vector<double>,std::vector<double>> ph_grad_val = Map2Physical(Jacobian, shp_fnc_dx, shp_fnc_dy);


											double bi_x = ph_grad_val.first[mycnt_1];
											double bi_y = ph_grad_val.second[mycnt_1];
											double bj_x = ph_grad_val.first[mycnt_2];
											double bj_y = ph_grad_val.second[mycnt_2];
											double wvol = gauss_x.second[g1] * gauss_y.second[g2];

											v += Jacobian.getDetVal() * (bi_x * bj_x + bi_y * bj_y) * wvol;
										}
									}
									int index1 = i1 * bspline_y->getNOF() + i2;
									int index2 = j1 * bspline_y->getNOF() + j2;
									double temp = A(index1, index2) + v;
									A.setValue(index1, index2, temp);

									mycnt_2++;
								}
							}
							mycnt_1++;
						}
					}
				}
			}
			else
			{
				if (trimmed_elements[cnt].second == 4) continue;
				int i_span_2 = bspline_y->findSpan(bspline_y->knots[ie2]);
				int mycnt_1 = 0;
				for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
				{
					for (int il_2 = 0; il_2 < bspline_y->getDegree() + 1; il_2++)
					{
						int mycnt_2 = 0;
						for (int jl_1 = 0; jl_1 < bspline_x->getDegree() + 1; jl_1++)
						{
							for (int jl_2 = 0; jl_2 < bspline_y->getDegree() + 1; jl_2++)
							{
								int i1 = i_span_1 - bspline_x->getDegree() + il_1;
								int j1 = i_span_1 - bspline_x->getDegree() + jl_1;

								int i2 = i_span_2 - bspline_y->getDegree() + il_2;
								int j2 = i_span_2 - bspline_y->getDegree() + jl_2;

								double v = 0.0;
								std::pair<std::vector<double>, std::vector<double>> gauss_x = bspline_x->calcGaussPts(bspline_x->getDegree() + 3, bspline_x->knots[ie1], bspline_x->knots[ie1 + 1]);
								std::pair<std::vector<double>, std::vector<double>> gauss_y = bspline_y->calcGaussPts(bspline_y->getDegree() + 3, bspline_y->knots[ie2], bspline_y->knots[ie2 + 1]);
								
								for (int g1 = 0; g1 < gauss_x.first.size(); g1++)
								{
									for (int g2 = 0; g2 < gauss_y.first.size(); g2++)
									{
										std::pair<std::vector<double>, std::vector<double>> eval_1 = bspline_x->eval(gauss_x.first[g1]);
										std::pair<std::vector<double>, std::vector<double>> eval_2 = bspline_y->eval(gauss_y.first[g2]);

										std::vector<double> shp_fnc_dx = createTensorProduct(eval_1.second, eval_2.first);
										std::vector<double> shp_fnc_dy = createTensorProduct(eval_1.first, eval_2.second);

										Matrix<double> Jacobian = calcJacobian(g1, g2, bspline_x->findSpan(gauss_x.first[g1]), bspline_y->findSpan(gauss_y.first[g2]), shp_fnc_dx, shp_fnc_dy);																		
										std::pair<std::vector<double>,std::vector<double>> ph_grad_val = Map2Physical(Jacobian, shp_fnc_dx, shp_fnc_dy);

										double bi_x = ph_grad_val.first[mycnt_1];
										double bi_y = ph_grad_val.second[mycnt_1];
										double bj_x = ph_grad_val.first[mycnt_2];
										double bj_y = ph_grad_val.second[mycnt_2];
										double wvol = gauss_x.second[g1] * gauss_y.second[g2];

										v += Jacobian.getDetVal() * (bi_x * bj_x + bi_y * bj_y) * wvol;																
									}
								}

								int index1 = i1 * bspline_y->getNOF() + i2;
								int index2 = j1 * bspline_y->getNOF() + j2;
								double temp = A(index1, index2) + v;
								A.setValue(index1, index2, temp);

								mycnt_2++;
							}
						}
						mycnt_1++;
					}
				}
			}
		}
	}
	stiff = A;
}

void Assembler_2D::calcRhs()
{
	// Assemble rhs vector
	std::vector<double> b(getNOF(), 0.0);
	int Nx = bspline_x->knots.size() - 1;
	int Ny = bspline_y->knots.size() - 1;
	int cnt = -1;
	for (int ie1 = 0; ie1 < Nx; ie1++)
	{
		int i_span_1 = bspline_x->findSpan(bspline_x->knots[ie1]);
		for (int ie2 = 0; ie2 < Ny; ie2++)
		{
			cnt++;
			if (trimmed_elements[cnt].first)
			{
				std::vector<std::vector<std::pair<double, double>>> triangles;
				if (trimmed_elements[cnt].second == 1)
				{
					triangles = construct_3_triangles(cnt);
				}
				else if (trimmed_elements[cnt].second == 2)
				{
					triangles = construct_2_triangles(cnt);
				}
				else if (trimmed_elements[cnt].second == 3)
				{
					triangles = construct_1_triangle(cnt);
				}
				for (auto triangle : triangles)
				{
					double a_x = std::min({triangle[0].first, triangle[1].first, triangle[2].first});
					double b_x = std::max({triangle[0].first, triangle[1].first, triangle[2].first});
					double a_y = std::min({triangle[0].second, triangle[1].second, triangle[2].second});
					double b_y = std::max({triangle[0].second, triangle[1].second, triangle[2].second});
					int i_span_2 = bspline_y->findSpan(bspline_y->knots[ie2]);
					int mycnt_1 = 0;
					for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
					{
						for (int il_2 = 0; il_2 < bspline_y->getDegree() + 1; il_2++)
						{
							int i1 = i_span_1 - bspline_x->getDegree() + il_1;
							int i2 = i_span_2 - bspline_y->getDegree() + il_2;

							double v = 0.0;
							std::pair<std::vector<double>, std::vector<double>> gauss_x = calcGaussPtsTria(a_x, b_x);
							std::pair<std::vector<double>, std::vector<double>> gauss_y = calcGaussPtsTria(a_y, b_y);
							for (int g1 = 0; g1 < gauss_x.first.size(); g1++)
							{
								for (int g2 = 0; g2 < gauss_y.first.size(); g2++)
								{
									if (g1 == (gauss_x.first.size()-1) && g2 == (gauss_y.first.size()-1)) break;
									std::pair<std::vector<double>, std::vector<double>> eval_1 = bspline_x->eval(gauss_x.first[g1]);
									std::pair<std::vector<double>, std::vector<double>> eval_2 = bspline_y->eval(gauss_y.first[g2]);

									std::vector<double> shp_fnc = createTensorProduct(eval_1.first, eval_2.first);
									std::vector<double> shp_fnc_dx = createTensorProduct(eval_1.second, eval_2.first);
									std::vector<double> shp_fnc_dy = createTensorProduct(eval_1.first, eval_2.second);

									Matrix<double> Jacobian = calcJacobian(g1, g2, bspline_x->findSpan(gauss_x.first[g1]), bspline_y->findSpan(gauss_y.first[g2]), shp_fnc_dx, shp_fnc_dy);

									std::pair<std::vector<double>, std::vector<double>> ph_val = Map2Physical(Jacobian, shp_fnc, shp_fnc);

									double bi_0 = ph_val.first[mycnt_1] * ph_val.second[mycnt_1];
									double wvol = gauss_x.second[g1] * gauss_y.second[g2];

									v += Jacobian.getDetVal() * bi_0 * f * wvol;
								}
							}
							int index = i1 * bspline_y->getNOF() + i2;
							b[index] += v;
						}
					}
					mycnt_1++;
				}
			}
			else
			{
				if (trimmed_elements[cnt].second == 4) continue;
				int i_span_2 = bspline_y->findSpan(bspline_y->knots[ie2]);
				int mycnt_1 = 0;

				for (int il_1 = 0; il_1 < bspline_x->getDegree() + 1; il_1++)
				{
					for (int il_2 = 0; il_2 < bspline_y->getDegree() + 1; il_2++)
					{
						int i1 = i_span_1 - bspline_x->getDegree() + il_1;
						int i2 = i_span_2 - bspline_y->getDegree() + il_2;

						double v = 0.0;
						std::pair<std::vector<double>, std::vector<double>> gauss_x = bspline_x->calcGaussPts(bspline_x->getDegree() + 3, bspline_x->knots[ie1], bspline_x->knots[ie1 + 1]);
						std::pair<std::vector<double>, std::vector<double>> gauss_y = bspline_y->calcGaussPts(bspline_y->getDegree() + 3, bspline_y->knots[ie2], bspline_y->knots[ie2 + 1]);
						for (int g1 = 0; g1 < gauss_x.first.size(); g1++)
						{
							for (int g2 = 0; g2 < gauss_y.first.size(); g2++)
							{
								std::pair<std::vector<double>, std::vector<double>> eval_1 = bspline_x->eval(gauss_x.first[g1]);
								std::pair<std::vector<double>, std::vector<double>> eval_2 = bspline_y->eval(gauss_y.first[g2]);

								std::vector<double> shp_fnc = createTensorProduct(eval_1.first, eval_2.first);
								std::vector<double> shp_fnc_dx = createTensorProduct(eval_1.second, eval_2.first);
								std::vector<double> shp_fnc_dy = createTensorProduct(eval_1.first, eval_2.second);

								Matrix<double> Jacobian = calcJacobian(g1, g2, bspline_x->findSpan(gauss_x.first[g1]), bspline_y->findSpan(gauss_y.first[g2]), shp_fnc_dx, shp_fnc_dy);	 			

								std::pair<std::vector<double>,std::vector<double>> ph_val = Map2Physical(Jacobian, shp_fnc, shp_fnc);

								double bi_0 = ph_val.first[mycnt_1] * ph_val.second[mycnt_1];
								double wvol = gauss_x.second[g1] * gauss_y.second[g2];

								v += Jacobian.getDetVal() * bi_0 * f * wvol;
							}
						}

						int index = i1 * bspline_y->getNOF() + i2;
						b[index] += v;
					}
				}
				mycnt_1++;
			}
		}
	}
	rhs = b;
}

std::pair<double, double> Assembler_2D::evaluate_trimming(double t)
{ 
	if (t < 0 || t >= 2.0*3.14159265) throw std::invalid_argument("Invalid parameter value");
	return std::make_pair(trimming[0] + trimming[2] * cos(t), trimming[1] + trimming[2] * sin(t));
}

std::pair<double, double> Assembler_2D::evaluate_trimming_der(double t)
{
	if (t < 0 || t >= 2.0*3.14159265) throw std::invalid_argument("Invalid parameter value");
	return std::make_pair(-trimming[2] * sin(t), trimming[2] * cos(t));
}

double Assembler_2D::projection_on_trimming(std::pair<double, double> point)
{
	std::vector<double> possible_projs;
	double u = 0.0;
	while (u <= 2 * 3.14159265)
	{
		std::pair<double, double> tangent_vec = evaluate_trimming_der(u);
		std::pair<double, double> di = std::make_pair(point.first - evaluate_trimming(u).first, point.second - evaluate_trimming(u).second);
		double dot = tangent_vec.first * di.first + tangent_vec.second * di.second;
		double d = sqrt(di.first * di.first + di.second * di.second);
		if (dot <= 1e-3 && dot >= -1e-3)
		{
			possible_projs.push_back(d);
		}
		u += 0.001;
	}

	double min = possible_projs[0];
	for (int i = 1; i < possible_projs.size(); i++)
	{
		if (possible_projs[i] < min) min = possible_projs[i];
	}

	return min;
}

int Assembler_2D::checkIfOutside(std::pair<double, double> point)
{
	std::vector<double> possible_projs;
	std::vector<std::pair<double, double>> possible_tangents;
	std::vector<std::pair<double, double>> possible_dis;
	double u = 0.0;
	while (u <= 2 * 3.14159265)
	{
		std::pair<double, double> tangent_vec = evaluate_trimming_der(u);
		std::pair<double, double> di = std::make_pair(point.first - evaluate_trimming(u).first, point.second - evaluate_trimming(u).second);
		double dot = tangent_vec.first * di.first + tangent_vec.second * di.second;
		double d = sqrt(di.first * di.first + di.second * di.second);
		if (dot <= 1e-3 && dot >= -1e-3)
		{
			possible_projs.push_back(d);
			possible_tangents.push_back(tangent_vec);
			possible_dis.push_back(di);
		}
		u += 0.001;
	}

	double min = 0;
	for (int i = 1; i < possible_projs.size(); i++)
	{
		if (possible_projs[i] < possible_projs[min])
			min = i;
	}

	std::vector<double> v1{possible_dis[min].first, possible_dis[min].second, 0.0};
	std::vector<double> v2{possible_tangents[min].first, possible_tangents[min].second, 0.0};
	double w1 = v1[1] * v2[2] - v1[2] * v2[1];
	double w2 = v1[2] * v2[0] - v1[0] * v2[2];
	double w3 = v1[0] * v2[1] - v1[1] * v2[0];

	if (w3 < 0 && w3 < -1e-5) return 1;
	else return 0;
}

void Assembler_2D::plot_trimming()
{
	std::string filename("trimming_curve.dat");
	std::ofstream my_file(filename);
	my_file << "variables= " << "\"x\"" << "," << "\"y\"" << "\n";
    my_file << "zone t= " << "\"1\"" << ",i=" << 63 << ",j=" << 63 << "\n";
	double u = 0.0;
	while (u <= 2 * 3.14159265)
	{
		std::pair<double, double> point = evaluate_trimming(u);
		my_file << point.first << " " << point.second << "\n";
		u += 0.1;
	}
	my_file.close();
}

std::vector<std::vector<std::pair<double, double>>> Assembler_2D::construct_3_triangles(int id)
{
	double min_x = std::min({elements_vertices[id][0].first, elements_vertices[id][1].first, elements_vertices[id][2].first, elements_vertices[id][3].first});
	double max_x = std::max({elements_vertices[id][0].first, elements_vertices[id][1].first, elements_vertices[id][2].first, elements_vertices[id][3].first});
	double min_y = std::min({elements_vertices[id][0].second, elements_vertices[id][1].second, elements_vertices[id][2].second, elements_vertices[id][3].second});
	double max_y = std::max({elements_vertices[id][0].second, elements_vertices[id][1].second, elements_vertices[id][2].second, elements_vertices[id][3].second});
	std::pair<double, double> diagonal, vertex3, vertex4;

	if ((trimmed_elements_info[id][0].first == trimmed_elements_info[id][1].first || trimmed_elements_info[id][0].first == trimmed_elements_info[id][2].first)
		&& (trimmed_elements_info[id][0].second == trimmed_elements_info[id][1].second || trimmed_elements_info[id][0].second == trimmed_elements_info[id][2].second))
		diagonal = trimmed_elements_info[id][0];
	else if ((trimmed_elements_info[id][1].first == trimmed_elements_info[id][0].first || trimmed_elements_info[id][1].first == trimmed_elements_info[id][2].first)
			&& (trimmed_elements_info[id][1].second == trimmed_elements_info[id][0].second || trimmed_elements_info[id][1].second == trimmed_elements_info[id][2].second))
			diagonal = trimmed_elements_info[id][1];
	else diagonal = trimmed_elements_info[id][2];

	for (int i = 0; i < 3; i++)
	{
		if (trimmed_elements_info[id][i].second == diagonal.second && trimmed_elements_info[id][i] != diagonal)
		{
			vertex3 = trimmed_elements_info[id][i];
		}
		if (trimmed_elements_info[id][i].first == diagonal.first && trimmed_elements_info[id][i] != diagonal)
		{
			vertex4 = trimmed_elements_info[id][i];
		}
	}

	std::pair<double, double> vertex1, vertex2;

	if (diagonal.first == min_x)
	{
		double s2 = max_x;
		double t = 0.0;
		double u = 0.0;
		while (u <= 2 * 3.14159265)
		{
			std::pair<double, double> point = evaluate_trimming(u);
			if (fabs(point.first - s2) < 1e-3)
			{
				if (point.second >= min_y && point.second <= max_y)
				{
					t = point.second;
					break;
				}
			}
			u += 0.001;
		}
		vertex1 = std::make_pair(s2, t);
		if (diagonal.second == min_y)
		{
			double t1 = max_y;
			double s = 0.0;
			u = 0.0;
			while (u <= 2 * 3.14159265)
			{
				std::pair<double, double> point = evaluate_trimming(u);
				if (fabs(point.second - t1) < 1e-3)
				{
					if (point.first >= min_x && point.first <= max_x)
					{
						s = point.first;
						break;
					}
				}
				u += 0.001;
			}
			vertex2 = std::make_pair(s, t1);
		}
		else if (diagonal.second == max_y)
		{
			double t1 = min_y;
			double s = 0.0;
			u = 0.0;
			while (u <= 2 * 3.14159265)
			{
				std::pair<double, double> point = evaluate_trimming(u);
				if (fabs(point.second - t1) < 1e-3)
				{
					if (point.first >= min_x && point.first <= max_x)
					{
						s = point.first;
						break;
					}
				}
				u += 0.001;
			}
			vertex2 = std::make_pair(s, t1);
		}
	}
	else if (diagonal.first == max_x)
	{
		double s2 = min_x;
		double t = 0.0;
		double u = 0.0;
		while (u <= 2 * 3.14159265)
		{
			std::pair<double, double> point = evaluate_trimming(u);
			if (fabs(point.first - s2) < 1e-3)
			{
				if (point.second >= min_y && point.second <= max_y)
				{
					t = point.second;
					break;
				}
			}
			u += 0.001;
		}
		vertex1 = std::make_pair(s2, t);
		if (diagonal.second == min_y)
		{
			double t1 = max_y;
			double s = 0.0;
			u = 0.0;
			while (u <= 2 * 3.14159265)
			{
				std::pair<double, double> point = evaluate_trimming(u);
				if (fabs(point.second - t1) < 1e-3)
				{
					if (point.first >= min_x && point.first <= max_x)
					{
						s = point.first;
						break;
					}
				}
				u += 0.001;
			}
			vertex2 = std::make_pair(s, t1);
		}
		else if (diagonal.second == max_y)
		{
			double t1 = min_y;
			double s = 0.0;
			u = 0.0;
			while (u <= 2 * 3.14159265)
			{
				std::pair<double, double> point = evaluate_trimming(u);
				if (fabs(point.second - t1) < 1e-3)
				{
					if (point.first >= min_x && point.first <= max_x)
					{
						s = point.first;
						break;
					}
				}
				u += 0.001;
			}
			vertex2 = std::make_pair(s, t1);
		}
	}

	std::vector<std::pair<double, double>> triangle1{vertex1, diagonal, vertex3};
	std::vector<std::pair<double, double>> triangle2{vertex1, diagonal, vertex2};
	std::vector<std::pair<double, double>> triangle3{vertex2, diagonal, vertex4};

	return std::vector<std::vector<std::pair<double, double>>>{triangle1, triangle2, triangle3};
}

std::vector<std::vector<std::pair<double, double>>> Assembler_2D::construct_2_triangles(int id)
{
	double min_x = std::min({elements_vertices[id][0].first, elements_vertices[id][1].first, elements_vertices[id][2].first, elements_vertices[id][3].first});
	double max_x = std::max({elements_vertices[id][0].first, elements_vertices[id][1].first, elements_vertices[id][2].first, elements_vertices[id][3].first});
	double min_y = std::min({elements_vertices[id][0].second, elements_vertices[id][1].second, elements_vertices[id][2].second, elements_vertices[id][3].second});
	double max_y = std::max({elements_vertices[id][0].second, elements_vertices[id][1].second, elements_vertices[id][2].second, elements_vertices[id][3].second});

	std::pair<double, double> vertex1, vertex2, vertex3;
	if (trimmed_elements_info[id][0].first == trimmed_elements_info[id][1].first)
	{
		if (trimmed_elements_info[id][0].second > trimmed_elements_info[id][1].second) vertex3 = trimmed_elements_info[id][0];
		else vertex3 = trimmed_elements_info[id][1];
		double t1 = max_y;
		double s = 0.0;
		double u = 0.0;
		while (u <= 2 * 3.14159265)
		{
			std::pair<double, double> point = evaluate_trimming(u);
			if (fabs(point.second - t1) < 1e-3)
			{
				if (point.first >= min_x && point.first <= max_x)
				{
					s = point.first;
					break;
				}
			}
			u += 0.001;
		}
		vertex1 = std::make_pair(s, t1);

		t1 = min_y;
		s = 0.0;
		u = 0.0;
		while (u <= 2 * 3.14159265)
		{
			std::pair<double, double> point = evaluate_trimming(u);
			if (fabs(point.second - t1) < 1e-3)
			{
				if (point.first >= min_x && point.first <= max_x)
				{
					s = point.first;
					break;
				}
			}
			u += 0.001;
		}
		vertex2 = std::make_pair(s, t1);
	}
	else if (trimmed_elements_info[id][0].second == trimmed_elements_info[id][1].second)
	{
		if (trimmed_elements_info[id][0].first > trimmed_elements_info[id][1].first) vertex3 = trimmed_elements_info[id][0];
		else vertex3 = trimmed_elements_info[id][1];
		double s1 = max_x;
		double t = 0.0;
		double u = 0.0;
		while (u <= 2 * 3.14159265)
		{
			std::pair<double, double> point = evaluate_trimming(u);
			if (fabs(point.first - s1) < 1e-3)
			{
				if (point.second >= min_y && point.second <= max_y)
				{
					t = point.second;
					break;
				}
			}
			u += 0.001;
		}
		vertex1 = std::make_pair(s1, t);

		s1 = min_x;
		t = 0.0;
		u = 0.0;
		while (u <= 2 * 3.14159265)
		{
			std::pair<double, double> point = evaluate_trimming(u);
			if (fabs(point.first - s1) < 1e-3)
			{
				if (point.second >= min_y && point.second <= max_y)
				{
					t = point.second;
					break;
				}
			}
			u += 0.001;
		}
		vertex2 = std::make_pair(s1, t);
	}

	std::vector<std::pair<double, double>> triangle1{vertex1, vertex3, vertex2};
	std::vector<std::pair<double, double>> triangle2{vertex2, trimmed_elements_info[id][0], trimmed_elements_info[id][1]};

	return std::vector<std::vector<std::pair<double, double>>>{triangle1, triangle2};
}

std::vector<std::vector<std::pair<double, double>>> Assembler_2D::construct_1_triangle(int id)
{
	double min_x = std::min({elements_vertices[id][0].first, elements_vertices[id][1].first, elements_vertices[id][2].first, elements_vertices[id][3].first});
	double max_x = std::max({elements_vertices[id][0].first, elements_vertices[id][1].first, elements_vertices[id][2].first, elements_vertices[id][3].first});
	double min_y = std::min({elements_vertices[id][0].second, elements_vertices[id][1].second, elements_vertices[id][2].second, elements_vertices[id][3].second});
	double max_y = std::max({elements_vertices[id][0].second, elements_vertices[id][1].second, elements_vertices[id][2].second, elements_vertices[id][3].second});
	double t1 = trimmed_elements_info[id][0].second;
	double s = 0.0;
	double u = 0.0;
	while (u <= 2 * 3.14159265)
	{
		std::pair<double, double> point = evaluate_trimming(u);
		if (fabs(point.second - t1) < 1e-3)
		{
			if (point.first >= min_x && point.first <= max_x)
			{
				s = point.first;
				break;
			}
		}
		u += 0.001;
	}
	std::pair<double, double> vertex1 = std::make_pair(s, t1);

	double s2 = trimmed_elements_info[id][0].first;
	double t = 0.0;
	u = 0.0;
	while (u <= 2 * 3.14159265)
	{
		std::pair<double, double> point = evaluate_trimming(u);
		if (fabs(point.first - s2) < 1e-3)
		{
			if (point.second >= min_y && point.second <= max_y)
			{
				t = point.second;
				break;
			}
		}
		u += 0.001;
	}
	std::pair<double, double> vertex2 = std::make_pair(s2, t);

	std::vector<std::pair<double, double>> triangle1{vertex1, trimmed_elements_info[id][0], vertex2};

	return std::vector<std::vector<std::pair<double, double>>>{triangle1};
}

std::pair<std::vector<double>, std::vector<double>> Assembler_2D::calcGaussPtsTria(double a, double b)
{
	std::vector<double> GS_pts_temp{1.0/6.0, 2.0/3.0};
	std::vector<double> GS_wgts{1.0/3.0, 1.0/3.0};
	// Convert to interval [a,b]
	std::vector<double> GS_pts;
	for (int i = 0; i < GS_pts_temp.size(); i++)
	{
		GS_pts.push_back(a * (1 - GS_pts_temp[i]) + b * GS_pts_temp[i]);
	}

	return std::make_pair(GS_pts, GS_wgts);
}
