#include "berendsen.h"
#include "statisticssampler.h"
#include <cmath>

BerendsenThermostat::BerendsenThermostat(double tbath, double relaxationTime, StatisticsSampler *statistics)
{
    m_tbath = tbath;
    m_relaxationTime = relaxationTime;
    m_statistics = statistics;
}

BerendsenThermostat::~BerendsenThermostat()
{
    delete m_statistics;
}

void BerendsenThermostat::scalingFactor(Atom *atom, double dt)
{
    atom->velocity * sqrt(1.0 + dt / m_relaxationTime * (m_tbath / m_statistics->temperature() - 1));
}
