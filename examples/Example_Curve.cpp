// Example of a NURBS curve (hemicycle)

#include <iostream>
#include "..\include\BsplineCurve.h"

int main()
{
    // - - - - - B-spline basis - - - - -
    int degree = 2;
    std::vector<double> valuesOfKnotVector{0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0};
    std::vector<double> weights{1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0};
    KnotVector<double> knotVector(degree, valuesOfKnotVector, weights);
    Bspline bspline_x(knotVector);

    // - - - - - B-spline curve - - - - -
    std::vector<Vertex<double>> controlPoints{{2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}};
    BsplineCurve curve(std::vector<Bspline>{bspline_x}, controlPoints);
    int resolution = 100;
    curve.plot(resolution, "curve.dat");

    return 0;
}