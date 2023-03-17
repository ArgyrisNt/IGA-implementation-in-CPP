#ifndef H_BOUNDCOND
#define H_BOUNDCOND

#include <iostream>
#include <vector>

class BoundCond
{
public:
    BoundCond() {}

    BoundCond(const std::pair<std::string, double> _west, const std::pair<std::string, double> _east,
              const std::pair<std::string, double> _north = std::make_pair("Dirichlet", 0.0), 
              const std::pair<std::string, double> _south = std::make_pair("Dirichlet", 0.0))
        : west(_west), east(_east), north(_north), south(_south) {}

    ~BoundCond() {}

    const std::pair<std::string, double> west;
    const std::pair<std::string, double> east;
    const std::pair<std::string, double> north;
    const std::pair<std::string, double> south;
};

#endif