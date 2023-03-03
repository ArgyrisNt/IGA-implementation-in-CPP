// Example of a NURBS curve

#include <iostream>
#include "..\IGA.h"

int main()
{
    // - - - - - B-spline basis - - - - -
    int p = 2;
    std::vector<double> U{0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0};
    std::vector<double> W{1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0};
    std::vector<std::vector<double>> ctrlPts{{2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}};
    Bspline bspline_x(p, U, W, ctrlPts);

    // - - - - - B-spline curve - - - - -
    BsplineCurve curve(bspline_x, ctrlPts);
    curve.evaluate();

    return 0;
}