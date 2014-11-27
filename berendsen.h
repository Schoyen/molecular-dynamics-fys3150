#pragma once
#include "atom.h"

class BerendsenThermostat
{
private:
    double m_tbath;
    double m_relaxationTime;

public:
    BerendsenThermostat(double tbath, double relaxationTime);
    ~BerendsenThermostat();
    void scalingFactor(Atom *atom, double temperature, double dt);
};
