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

double BerendsenThermostat::scalingFactor(double temperature, double dt)
{
    return sqrt(1.0 + dt / m_relaxationTime * (m_tbath / temperature - 1));
}
