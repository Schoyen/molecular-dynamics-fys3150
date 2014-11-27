#include "potential.h"

Potential::Potential() :
    m_potentialEnergy(0)
{

}

double Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

double Potential::pressure()
{
    return m_pressure;
}
