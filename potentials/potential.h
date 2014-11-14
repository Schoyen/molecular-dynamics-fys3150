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
    BerendsenThermostat *m_berendsen;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    void addThermostat(BerendsenThermostat *berendsen);
    double potentialEnergy();
    double kineticEnergy();
    double temperature();
    BerendsenThermostat *berendsen() {return m_berendsen;}
};
