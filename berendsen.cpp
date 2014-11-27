#include "berendsen.h"
#include <cmath>

BerendsenThermostat::BerendsenThermostat(double tbath, double relaxationTime)
{
    m_tbath = tbath;
    m_relaxationTime = relaxationTime;
}

BerendsenThermostat::~BerendsenThermostat()
{

}

void BerendsenThermostat::scalingFactor(Atom *atom, double temperature, double dt)
{
    atom->velocity * sqrt(1.0 + dt / m_relaxationTime * (m_tbath / temperature - 1));
}
