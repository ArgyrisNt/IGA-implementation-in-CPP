#include "..\include\MultiBspline.h"
#include <iostream>
#include <string>

void MultiBspline::plot(int resolution, const std::string &filename)
{
    for (int i = 0; i < dimension; ++i)
        bsplines[i].plot(resolution, std::to_string(i + 1) + filename);
}