#pragma once

#include "..\include\KnotVector.h"
#include <iostream>
#include <vector>

class Bspline
{
public:
	Bspline(int newDegree, KnotVector<double>& knewKnotVector, std::vector<double> &newWeights);
	Bspline(const Bspline&);

	~Bspline();

	Bspline& operator=(const Bspline&);

	std::pair<std::vector<double>, std::vector<double>> GaussPointsAndWeights(int, const double, const double); 
	std::pair<std::vector<double>, std::vector<double>> evaluateAtPoint(const double value, bool all = false); // all basis functions if all = true
	void plot(int resolution);

	void setWeights(std::vector<double> new_weights);
	void setKnotvector(KnotVector<double>& newKnotVector) { knotVector = newKnotVector; }

	KnotVector<double>& getKnotvector() { return knotVector; };
	const int getDegree() const { return degree; }
	const int getNumberOfBasisFunctions() const { return numberOfBasisFunctions; }
	std::vector<double>& getWeights() { return weights; }

private:
	void basisFunctionsOfDegree(int level, double value, std::vector<double> &values, std::vector<double>& derivatives);
	void computeActiveBasisFunctions(double value, std::vector<double> &values, std::vector<double> &derivatives);
	std::vector<double> initializeBasisFunctions(double value);
	void mapValuesToDomain(std::vector<double> &GaussPoints, const double left, const double right);

	KnotVector<double> knotVector;
	std::vector<double> weights;
	int degree;
	int numberOfBasisFunctions;
};

class BsplineCurve
{
public:
	BsplineCurve(Bspline new_bspline_x, std::vector<std::vector<double>> new_controlPoints) : bspline_x(new_bspline_x), controlPoints(new_controlPoints) {}

	~BsplineCurve() {}

	std::pair<double, double> evaluateAtPoint(double point);
	void plot(int resolution);

	Bspline &getBspline_x() { return bspline_x; }
	std::vector<std::vector<double>> &getControlPoints() { return controlPoints; };

private:
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

	std::pair<double, double> evaluateAtPoint(std::pair<double, double>&& point);
	void plot(int resolution);
	void uniformRefine_x();
	void uniformRefine_y();
	void knotInsertion(KnotVector<double> &vector, std::vector<std::vector<double>> &points, double newKnot);

private:
	std::vector<std::vector<double>> pointsOfParametricCurve(int direction, int level); // direction = 0 for x
	void refineParametricCurve(KnotVector<double> &vector, std::vector<std::vector<double>> &points);

	Bspline bspline_x;
	Bspline bspline_y;
	std::vector<std::vector<double>> controlPoints;
};