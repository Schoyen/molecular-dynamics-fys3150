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

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    double x, y, z;
    double lengthOfCube = 5;
    for (int n = 0; n < 100; n++) {
        x = m_atoms[n]->position.x();
        y = m_atoms[n]->position.y();
        z = m_atoms[n]->position.z();
        if (x < -lengthOfCube) {
            x = m_atoms[n]->position.x() + lengthOfCube;
        }
        else if (x >= lengthOfCube) {
            x = m_atoms[n]->position.x() - lengthOfCube;
        }
        if (y < -lengthOfCube) {
            y = m_atoms[n]->position.y() + lengthOfCube;
        }
        else if (y >= lengthOfCube) {
            y = m_atoms[n]->position.y() - lengthOfCube;
        }
        if (z < -lengthOfCube) {
            z = m_atoms[n]->position.z() + lengthOfCube;
        }
        else if (z >= lengthOfCube) {
            z = m_atoms[n]->position.z() - lengthOfCube;
        }
        m_atoms[n]->position = vec3(x, y, z);
    }
}

void System::removeMomentum() {
    // Initially, when the atoms are given random velocities, there is a non-zero net momentum. We don't want any drift in the system, so we need to remove it.
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
