#pragma once
#include "integrator.h"

class VelocityVerlet : public Integrator
{
private:
    void firstKick(System *system, double dt);
    void halfKick(System *system, double dt);
    void move(System *system, double dt);
    bool m_firstStep;
    bool m_thermostat;
public:
    VelocityVerlet();
    ~VelocityVerlet();
    virtual void integrate(System *system, double dt, bool thermostatOn);
};
