#pragma once

#include <iostream>
#include <vector>

class Bspline
{
public:
	// Constructors
	Bspline(); // default constructor
	Bspline(int, std::vector<double>&, std::vector<double>&); // initialize a basis by defining a degree and a knotvector
	Bspline(const Bspline&); // copy constructor
	Bspline(const double start, const double end, int new_degree, int numberOfElements, std::vector<double>& new_weights); // initialize a basis by defining the start and the end of the knotvector, the degree and the number of elements

	// Destructor
	~Bspline();

	// Overload operator
	Bspline& operator=(const Bspline&);

	// Member functions
	std::pair<std::vector<double>, std::vector<double>> GaussPointsAndWeights(int, const double, const double); // compute Gauss points and their weigths
	std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(const double value, bool all = false); // eval all active basis functions at x (all basis functions if all=true)
	void plot(int resolution);
	int findSpanInVector(const double value, std::vector<double> vector = {});

	// Member setter functions
	void setWeights(std::vector<double> new_weights);
	void setKnotvector(std::vector<double> new_knotvector) { knotvector = new_knotvector; }

	// Member getter functions
	std::vector<double>& getKnotvector() { return knotvector; };
	const int getDegree() const { return degree; }
	const int getNumberOfBasisFunctions() const { return numberOfBasisFunctions; }
	std::vector<double>& getWeights() { return weights; }

	// Member variables
	std::vector<double> distinctKnots;

private:
	// Member local functions
	void computeDistinctKnots(); // compute discrete knots
	void basisFunctionsOfDegree(int level, double value, std::vector<double> &values, std::vector<double>& derivatives);
	void computeActiveBasisFunctions(double value, std::vector<double> &values, std::vector<double> &derivatives);
	std::vector<double> initializeBasisFunctions(double value);

	// Member variables
	std::vector<double> knotvector;
	std::vector<double> weights;
	int degree;
	int numberOfBasisFunctions;
};

class BsplineCurve
{
public:
	BsplineCurve(Bspline new_bspline_x, std::vector<std::vector<double>> new_controlPoints) : bspline_x(new_bspline_x), controlPoints(new_controlPoints) {}

	~BsplineCurve() {}

	void plot(int resolution);

	Bspline bspline_x;
	std::vector<std::vector<double>> controlPoints;
};

class BsplineSurface
{
public:
	BsplineSurface(const Bspline &new_bspline_x, const Bspline &new_bspline_y, std::vector<std::vector<double>> new_controlPoints)
		: bspline_x(new_bspline_x), bspline_y(new_bspline_y), controlPoints(new_controlPoints) {}

	~BsplineSurface() {}

	BsplineSurface &operator=(const BsplineSurface &);

	Bspline &getBspline_x() { return bspline_x; }
	Bspline &getBspline_y() { return bspline_y; }
	std::vector<std::vector<double>> &getControlPoints() { return controlPoints; };

	void setControlPoints(std::vector<std::vector<double>> &new_controlPoints) { controlPoints = new_controlPoints; }
	void setBspline_x(Bspline &new_bspline_x) { bspline_x = new_bspline_x; }
	void setBspline_y(Bspline &new_bspline_y) { bspline_y = new_bspline_y; }

	void plot(int resolution);
	void uniformRefine_x();
	void uniformRefine_y();
	void knotInsertion(std::vector<double> &vector, std::vector<std::vector<double>> &points, double newKnot);

private:
	void refineParametricCurve(std::vector<double> &vector, std::vector<std::vector<double>> &points);
	std::vector<std::vector<double>> pointsOfParametricCurve(int direction, int level); // direction = 0 for x

	Bspline bspline_x;
	Bspline bspline_y;
	std::vector<std::vector<double>> controlPoints;
};