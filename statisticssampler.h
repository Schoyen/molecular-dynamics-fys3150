#pragma once
#include "system.h"
#include "unitconverter.h"
#include <fstream>

using namespace std;

class StatisticsSampler
{
private:
    double m_netMomentum;
    double m_netMomentumAfter;
    double m_kineticEnergy;
    double m_potentialEnergy;
    double m_temperature;
    double m_numberDensity;
    double m_pressure;
    ofstream m_kineticEnergyFile;
    ofstream m_potentialEnergyFile;
    ofstream m_netMomentumFile;
    ofstream m_temperatureFile;
    ofstream m_numberDensityFile;
    ofstream m_pressureFile;
    ofstream m_heatCapacityFile;
    // Add optional pressure minus heat capacity.
    bool m_sampledNetMomentum;
public:
    StatisticsSampler();
    ~StatisticsSampler();
    void createFiles();
    void closeFiles();
    void sample(System *system, int timestep);
    void sampleKineticEnergy(System *system);
    void samplePotentialEnergy(System *system);
    void sampleNetMomentum(System *system);
    void sampleTemperature(System *system);
    void sampleNumberDensity(System *system);
    void samplePressure(System *system);
    double totalEnergy() {return m_kineticEnergy + m_potentialEnergy;}
};
