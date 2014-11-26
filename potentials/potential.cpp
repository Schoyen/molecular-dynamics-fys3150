#include "potential.h"

Potential::Potential() :
    m_potentialEnergy(0),
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

double Potential::pressure()
{
    return m_pressure;
}
