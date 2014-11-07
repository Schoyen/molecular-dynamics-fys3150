#pragma once
#include "system.h"
#include "unitconverter.h"
#include <iostream>

class StatisticsSampler
{
private:
    double m_netMomentum;
    double m_netMomentumAfter;
    double m_kineticEnergy;
    double m_potentialEnergy;
    double m_temperature;
    bool sampledNetMomentum;
public:
    StatisticsSampler();
    ~StatisticsSampler();
    void sample(System *system, std::string filename);
    void sampleKineticEnergy(System *system);
    void samplePotentialEnergy(System *system);
    void sampleNetMomentum(System *system);
    void sampleTemperature(System *system);
};
