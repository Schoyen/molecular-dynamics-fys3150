#pragma once
#include "statisticssampler.h"

class Atom;

class BerendsenThermostat
{
private:
    double m_tbath;
    double m_relaxationTime;
    StatisticsSampler *m_statistics;

public:
    BerendsenThermostat(double tbath, double relaxationTime, StatisticsSampler *statistics);
    ~BerendsenThermostat();
    void scalingFactor(Atom *atom, double dt);
};
