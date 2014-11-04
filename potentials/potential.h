#pragma once
#include <string>
#include <vector>
#include "../system.h"

class Potential
{
protected:
    double m_potentialEnergy;
    double m_kineticEnergy; // Added for use in the calculateForces-function.
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    double potentialEnergy();
    void setPotentialEnergy(double potentialEnergy);
};
