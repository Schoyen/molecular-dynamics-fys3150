#pragma once
#include <string>
#include <vector>
#include "../system.h"

class System; class BerendsenThermostat;

class Potential
{
protected:
    double m_potentialEnergy;
    double m_kineticEnergy;
    double m_temperature;
    double m_numberDensity;
    double m_pressure;
    BerendsenThermostat *m_berendsen;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    void addThermostat(BerendsenThermostat *berendsen);
    double potentialEnergy();
    double kineticEnergy();
    double temperature();
    double numberDensity();
    double pressure();
    BerendsenThermostat *berendsen() {return m_berendsen;}
};
