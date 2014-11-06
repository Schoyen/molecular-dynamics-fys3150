#pragma once
#include "system.h"

/*
 * Add samples.
 */
class StatisticsSampler
{
public:
    StatisticsSampler();
    ~StatisticsSampler();
    void sample(System *system);
    void sampleKineticEnergy(System *system);
};
