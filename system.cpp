#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "unitconverter.h"

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0)
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    m_atoms.clear();
}

/*
 * The atoms are not properly contained inside the box. Is this correct?
 * Do I need to convert to another system?
 */
void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    double x, y, z;
    for (int n = 0; n < 100; n++) {
        x = m_atoms[n]->position.x();
        y = m_atoms[n]->position.y();
        z = m_atoms[n]->position.z();
        if (x < -m_systemSize.x() * 0.5) {
            x += (m_systemSize.x() * 0.5);
        }
        else if (x >= m_systemSize.x() * 0.5) {
            x -= (m_systemSize.x() * 0.5);
        }
        if (y < -m_systemSize.y() * 0.5) {
            y += (m_systemSize.y() * 0.5);
        }
        else if (y >= m_systemSize.y() * 0.5) {
            y -= (m_systemSize.y() * 0.5);
        }
        if (z < -m_systemSize.z() * 0.5) {
            z += (m_systemSize.z() * 0.5);
        }
        else if (z >= m_systemSize.z() * 0.5) {
            z -= (m_systemSize.z() * 0.5);
        }
        m_atoms[n]->position = vec3(x, y, z);
    }
}

/*
 * Is this enough?
 * Do I need to convert to another system?
 */
void System::removeMomentum() {
    // Initially, when the atoms are given random velocities, there is a non-zero net momentum. We don't want any drift in the system, so we need to remove it.
    
    // First we calculate the net momentum of the system.
    // We then subtract this sum from the velocity of each atom.
    
    vec3 netMomentum (0.0, 0.0, 0.0);
    vec3 speedOfCM (0.0, 0.0, 0.0);
    vec3 temp;
    double totalMass = 0;
    for (int n = 0; n < 100; n++) {
        temp = m_atoms[n]->velocity*(m_atoms[n]->mass());
        speedOfCM = speedOfCM + temp;
        totalMass += m_atoms[n]->mass();
    }

    speedOfCM = speedOfCM/totalMass;
    for (int n = 0; n < 100; n++) {
        m_atoms[n]->velocity = m_atoms[n]->velocity - speedOfCM;
    }
}

void System::resetForcesOnAllAtoms() {

}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant) {

}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}
