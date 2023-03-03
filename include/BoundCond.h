#pragma once

#include <iostream>
#include <vector>

class BoundCond
{
public:
    // Constructors
    BoundCond() {}
    BoundCond(const std::string _west, const std::string _east, double _wval, double _eval)
        : west(_west), east(_east), wval(_wval), eval(_eval) {};
        
    BoundCond(const std::string _west, const std::string _east, const std::string _north, const std::string _south, double _wval, double _eval, double _nval, double _sval)
        : west(_west), east(_east), north(_north), south(_south), wval(_wval), eval(_eval), nval(_nval), sval(_sval) {};

    // Destructor
    ~BoundCond() {}

    // Member getter functions
    const std::string& getWest() { return west; }
    const std::string& getEast() { return east; }
    const std::string& getNorth() { return north; }
    const std::string& getSouth() { return south; }
    double getWval() { return wval; }
    double getEval() { return eval; }
    double getNval() { return nval; }
    double getSval() { return sval; }

private:
    // Member variables
    const std::string west;
    const std::string east;
    const std::string north;
    const std::string south;
    double wval;
    double eval;
    double nval;
    double sval;
};