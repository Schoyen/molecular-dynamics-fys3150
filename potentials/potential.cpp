#include "potential.h"

Potential::Potential() :
    m_potentialEnergy(0)
{

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
