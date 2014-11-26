#pragma once
#include <string>
#include <vector>
#include "../system.h"

class System; class BerendsenThermostat;

class Potential
{
protected:
    double m_potentialEnergy;
    double m_pressure; // Sum of all forces.
    BerendsenThermostat *m_berendsen;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    virtual void calculateForcesOld(System *system) = 0;
    void addThermostat(BerendsenThermostat *berendsen);
    double potentialEnergy();
    double pressure();
    BerendsenThermostat *berendsen() {return m_berendsen;}
};
