// Example of a NURBS curve

#include <iostream>
#include "..\IGA.h"

int main()
{
    // - - - - - B-spline basis - - - - -
    int degree = 2;
    std::vector<double> knotVector{0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0};
    std::vector<double> weights{1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0};
    Bspline bspline_x(degree, knotVector, weights);

    // - - - - - B-spline curve - - - - -
    std::vector<std::vector<double>> controlPoints{{2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}};
    BsplineCurve curve(bspline_x, controlPoints);
    int resolution = 100;
    curve.plot(resolution);

    return 0;
}