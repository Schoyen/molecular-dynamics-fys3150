#include "velocityverlet.h"
#include "../system.h"
#include "../berendsen.h"

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true) // This will set the variable m_firstStep to false when the object is created
{

}

VelocityVerlet::~VelocityVerlet()
{

}

void VelocityVerlet::firstKick(System *system, double dt)
{
    m_firstStep = false;
    system->calculateForces();
    halfKick(system, dt);
}

void VelocityVerlet::halfKick(System *system, double dt)
{
    if (m_thermostat) {
        for (int n = 0; n < (int) system->atoms().size(); n++) {
            Atom *atom = system->atoms()[n];
            double timestepDividedByTwoTimesMass = 0.5 * dt / atom->mass();
            atom->velocity.addAndMultiply(atom->force, timestepDividedByTwoTimesMass); // v += F/(2*m)*dt
            system->berendsen()->scalingFactor(atom, system->temperature, dt);
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

    // This does not seem to speed things up...
    std::cout << "YAY2" << std::endl;
    for (int n = 0; n < (int) system->atoms().size(); n++) {
        Atom *atom = system->atoms()[n];
        atom->position.addAndMultiply(atom->velocity, dt);

        // Idiotic boundary conditions.
        if (atom->position.x() < 0) atom->position[0] += system->systemSize().x();
        else if (atom->position.x() >= system->systemSize().x()) atom->position[0] -= system->systemSize().x();
        if (atom->position.y() < 0) atom->position[1] += system->systemSize().y();
        else if (atom->position.y() >= system->systemSize().y()) atom->position[1] -= system->systemSize().y();
        if (atom->position.z() < 0) atom->position[2] += system->systemSize().z();
        else if (atom->position.z() >= system->systemSize().z()) atom->position[2] -= system->systemSize().z();

        // Putting atoms in correct cells.
        int cx = int(atom->position.x() / system->systemSize().x() * system->numberOfCellsX);
        int cy = int(atom->position.y() / system->systemSize().y() * system->numberOfCellsY);
        int cz = int(atom->position.z() / system->systemSize().z() * system->numberOfCellsZ);
        system->celllist()->listOfCells()[cx * system->numberOfCellsY * system->numberOfCellsZ + cy * system->numberOfCellsZ + cz]->addAtom(atom);
    }
    std::cout << "YAY2" << std::endl;
}

void VelocityVerlet::integrate(System *system, double dt, bool thermostatOn)
{
    m_thermostat = thermostatOn;
    if(m_firstStep) firstKick(system, dt);
    else halfKick(system, dt);
    move(system, dt);
    // Attempting to use periodic boundary conditions and the calculation of cell position of atoms inside the move method.
    //system->applyPeriodicBoundaryConditions();
    system->calculateForces();
    halfKick(system, dt);
}
