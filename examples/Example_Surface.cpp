// Example of a NURBS surface

#include <iostream>
#include "..\include\BsplineSurface.h"

int main()
{
    // - - - - - B-spline basis on x-direction - - - - - 
    int x_degree = 2;
    std::vector<double> x_values{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    std::vector<double> x_weights{1.0, 1.0, 1.0};
    KnotVector<double> x_knotVector(x_degree, x_values, x_weights);
    Bspline bspline_x(x_knotVector);

    // - - - - - B-spline basis on y-direction - - - - - 
    int y_degree = 2;
    std::vector<double> y_values{0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    std::vector<double> y_weights{1.0, 1.0, 1.0, 1.0};
    KnotVector<double> y_knotVector(y_degree, y_values, y_weights);
    Bspline bspline_y(y_knotVector);

    // - - - - - Assempler info - - - - -
    std::vector<std::vector<double>> controlPoints{ {0.0, 0.0}, {0.0, 4.0}, {4.0, 8.0}, {8.0, 8.0},
                                              {2.0, 0.0}, {2.0, 3.0}, {5.0, 6.0}, {8.0, 6.0},
                                              {4.0, 0.0}, {4.0, 2.0}, {6.0, 4.0}, {8.0, 4.0} };

    // - - - - - B-spline surface - - - - -
    BsplineSurface surface(bspline_x, bspline_y, controlPoints);

    int numberOfRefinements = 3;
    for (int i = 0; i < numberOfRefinements; i++)
    {
        surface.uniformRefine_x();
        surface.uniformRefine_y();
    }

    int resolution = 100;
    surface.plot(resolution);

    return 0;
}