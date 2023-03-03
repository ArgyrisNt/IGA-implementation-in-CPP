// Example of a NURBS surface

#include <iostream>
#include "..\IGA.h"

int main()
{
    // - - - - - B-spline basis on x-direction - - - - - 
    int p = 2;
    std::vector<double> U{ 0.0,0.0,0.0,1.0,1.0,1.0 };
    std::vector<double> W{ 1.0,1.0,1.0 };
    Bspline bspline_x(p, U, W);

    // - - - - - B-spline basis on y-direction - - - - - 
    int q = 2;
    std::vector<double> V{ 0.0,0.0,0.0,0.5,1.0,1.0,1.0 };
    std::vector<double> Q{ 1.0,1.0,1.0,1.0 };
    Bspline bspline_y(q, V, Q);

    // - - - - - Assempler info - - - - -
    std::vector<std::vector<double>> ctrlPts{ {0.0, 0.0}, {0.0, 4.0}, {4.0, 8.0}, {8.0, 8.0},
                                              {2.0, 0.0}, {2.0, 3.0}, {5.0, 6.0}, {8.0, 6.0},
                                              {4.0, 0.0}, {4.0, 2.0}, {6.0, 4.0}, {8.0, 4.0} };

    // - - - - - B-spline surface - - - - -
    BsplineSurface surface(bspline_x, bspline_y, ctrlPts);
    std::cout << surface.getCtrlPts().size() << std::endl;
    surface.uniformRefine_x();
    surface.evaluate();

    return 0;
}