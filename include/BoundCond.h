#pragma once

#include <iostream>
#include <vector>

class BoundCond
{
public:
    // Constructors
    BoundCond() {}
    BoundCond(const std::string newWestType, const std::string newEastType, double newWestValue, double newEastValue)
        : westType(newWestType), eastType(newEastType), westValue(newWestValue), eastvalue(newEastValue){};

    BoundCond(const std::string newWestType, const std::string newEastType, const std::string newNorthType, const std::string newSouthType, double newWestValue, double newEastValue, double newNorthValue, double newSouthValue)
        : westType(newWestType), eastType(newEastType), northType(newNorthType), southType(newSouthType), westValue(newWestValue), eastvalue(newEastValue), northValue(newNorthValue), southValue(newSouthValue){};

    // Destructor
    ~BoundCond() {}

    // Member getter functions
    const std::string &getWestType() { return westType; }
    const std::string &getEastType() { return eastType; }
    const std::string &getNorthType() { return northType; }
    const std::string &getSouthType() { return southType; }
    double getWestValue() { return westValue; }
    double getEastValue() { return eastvalue; }
    double getNorthValue() { return northValue; }
    double getSouthValue() { return southValue; }

private:
    // Member variables
    const std::string westType;
    const std::string eastType;
    const std::string northType;
    const std::string southType;
    double westValue;
    double eastvalue;
    double northValue;
    double southValue;
};