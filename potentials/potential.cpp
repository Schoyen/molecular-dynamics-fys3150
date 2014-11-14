#include "potential.h"

Potential::Potential() :
    m_potentialEnergy(0),
    m_kineticEnergy(0),
    m_temperature(0),
    m_berendsen(0)
{

}

void Potential::addThermostat(BerendsenThermostat *berendsen)
{
    m_berendsen = berendsen;
}

double Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

double Potential::kineticEnergy()
{
    return m_kineticEnergy;
}

double Potential::temperature()
{
    return m_temperature;
}
