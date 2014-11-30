#pragma once
#include "unitconverter.h"
#include "system.h"
#include <fstream>

using namespace std;

class StatisticsSampler
{
private:
    double m_netMomentum;
    double m_netMomentumAfter;
    double m_kineticEnergy;
    double m_potentialEnergy;
    double m_totalEnergy;
    double m_temperature;
    double m_numberDensity;
    double m_pressure;
    double m_time;
    double m_heatCapacity;
    ofstream m_kineticEnergyFile;
    ofstream m_potentialEnergyFile;
    ofstream m_totalEnergyFile;
    ofstream m_netMomentumFile;
    ofstream m_temperatureFile;
    ofstream m_numberDensityFile;
    ofstream m_pressureFile;
    ofstream m_heatCapacityFile;
    ofstream m_timeFile;
    // Add optional pressure minus heat capacity.
    bool m_sampledNetMomentum;
    bool m_sampledNumberDensity;
public:
    static double m_kineticEnergySquared;
    static double m_totalKineticEnergy;
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
    void sampleTotalEnergy(System *system);
    void sampleTime(System *system);
    void sampleKineticEnergySquared(System *system);
    void sampleTotalKineticEnergy(System *system);
    void sampleHeatCapacity(System *system);
    double temperature() {return m_temperature;}
};
