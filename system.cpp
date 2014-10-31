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
    double xlen = m_systemSize.x();
    double ylen = m_systemSize.y();
    double zlen = m_systemSize.z();
    for (int n = 0; n < (int) m_atoms.size(); n++) {
        x = m_atoms[n]->position.x();
        y = m_atoms[n]->position.y();
        z = m_atoms[n]->position.z();
        if (x < -xlen * 0.5) {
            x += xlen;
        }
        else if (x >= xlen * 0.5) {
            x -= xlen;
        }
        if (y < -ylen * 0.5) {
            y += ylen;
        }
        else if (y >= ylen * 0.5) {
            y -= ylen;
        }
        if (z < -zlen * 0.5) {
            z += zlen;
        }
        else if (z >= zlen * 0.5) {
            z -= zlen;
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
    
    vec3 speedOfCM (0.0, 0.0, 0.0);
    vec3 temp;
    double totalMass = 0;
    for (int n = 0; n < (int) m_atoms.size(); n++) {
        temp = m_atoms[n]->velocity*(m_atoms[n]->mass());
        speedOfCM = speedOfCM + temp;
        totalMass += m_atoms[n]->mass();
    }

    speedOfCM = speedOfCM/totalMass;
    for (int n = 0; n < (int) m_atoms.size(); n++) {
        m_atoms[n]->velocity = m_atoms[n]->velocity - speedOfCM;
    }

}

void System::resetForcesOnAllAtoms() {

}

/*
 * Is this correct?
 */
void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant) {
    int totalNumberOfUnitCells = numberOfUnitCellsEachDimension * numberOfUnitCellsEachDimension * numberOfUnitCellsEachDimension;
    //int numberOfAtomsPerUnitCell = 4; // For Argon.
    vec3 r1 = vec3(0.0, 0.0, 0.0);
    vec3 r2 = vec3(latticeConstant/2.0, latticeConstant/2.0, 0.0);
    vec3 r3 = vec3(0, latticeConstant/2.0, latticeConstant/2.0);
    vec3 r4 = vec3(latticeConstant/2.0, 0.0, latticeConstant/2.0);
    //vec3 UNIT = vec3(latticeConstant, latticeConstant, latticeConstant);
    vector<vec3> R;
    vec3 temp;
    // Creating positions for the center of each unit cell.
    for (int i = 0; i < numberOfUnitCellsEachDimension; i++) {
        for (int j = 0; j < numberOfUnitCellsEachDimension; j++) {
            for (int k = 0; k < numberOfUnitCellsEachDimension; k++) {
                temp = vec3(i*latticeConstant, j*latticeConstant, k*latticeConstant);
                R.push_back(temp);
            }
        }
    }

    for (int n = 0; n < totalNumberOfUnitCells; n++) {
        Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom1->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
        atom1->position = r1 + R[n];
        m_atoms.push_back(atom1);

        Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom2->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
        atom2->position = r2 + R[n];
        m_atoms.push_back(atom2);

        Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom3->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
        atom3->position = r3 + R[n];
        m_atoms.push_back(atom3);

        Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom4->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
        atom4->position = r4 + R[n];
        m_atoms.push_back(atom4);
    }
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
