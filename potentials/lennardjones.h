#pragma once
#include "potential.h"

class LennardJones : public Potential
{
private:
    double m_sigma;
    double m_epsilon;
public:
    LennardJones(double sigma, double epsilon, BerendsenThermostat *berendsen);
    ~LennardJones() {}
    virtual void calculateForces(System *system);
};
