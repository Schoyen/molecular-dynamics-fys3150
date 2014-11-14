#include "velocityverlet.h"
#include "../system.h"

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true) // This will set the variable m_firstStep to false when the object is created
{

}

VelocityVerlet::~VelocityVerlet()
{

}

void VelocityVerlet::firstKick(System *system, double dt, bool thermostatOn)
{
    m_firstStep = false;
    system->calculateForces();
    halfKick(system, dt, thermostatOn);
}

void VelocityVerlet::halfKick(System *system, double dt, bool thermostatOn)
{
    if (thermostatOn) {
        for (int n = 0; n < (int) system->atoms().size(); n++) {
            Atom *atom = system->atoms()[n];
            double timestepDividedByTwoTimesMass = dt / (2.0*atom->mass());
            atom->velocity.addAndMultiply(atom->force, timestepDividedByTwoTimesMass); // v += F/(2*m)*dt
            atom->velocity = atom->velocity * system->potential()->berendsen()->scalingFactor(system->potential()->temperature(), dt);
        }
    } else {
        for (int n = 0; n < (int) system->atoms().size(); n++) {
            Atom *atom = system->atoms()[n];
            double timestepDividedByTwoTimesMass = dt / (2.0*atom->mass());
            atom->velocity.addAndMultiply(atom->force, timestepDividedByTwoTimesMass); // v += F/(2*m)*dt
        }
    }
}

void VelocityVerlet::move(System *system, double dt)
{
    for (int n = 0; n < (int) system->atoms().size(); n++) {
        Atom *atom = system->atoms()[n];
        atom->position.addAndMultiply(atom->velocity, dt);
    }
}

void VelocityVerlet::integrate(System *system, double dt, bool thermostatOn)
{
    if (thermostatOn) {
        if(m_firstStep) firstKick(system, dt);
        else halfKick(system, dt);
        move(system, dt);
        system->applyPeriodicBoundaryConditions();
        system->calculateForces();
        halfKick(system, dt);
        system->resetForcesOnAllAtoms();
    } else {
        if(m_firstStep) firstKick(system, dt);
        else halfKick(system, dt);
        move(system, dt);
        system->applyPeriodicBoundaryConditions();
        system->calculateForces();
        halfKick(system, dt);
        system->resetForcesOnAllAtoms();
    }
}
