#pragma once
#include <string>
#include <vector>
#include "../system.h"

class System;

class Potential
{
protected:
    double m_potentialEnergy;
    double m_pressure; // Sum of all forces.
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    virtual void calculateForcesOld(System *system) = 0;
    double potentialEnergy();
    double pressure();
};
