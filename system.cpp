#include "system.h"
#include "integrators/integrator.h"
#include "unitconverter.h"
#include "berendsen.h"
#include <fstream>

using namespace std;

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_rcut(0),
    m_steps(0),
    //m_thermostat(0),
    m_celllist(0)
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    delete m_celllist;
    //delete m_thermostat;
    m_atoms.clear();
}

void System::minimumImageCriterion(vec3 &pos)
{
    if (pos.x() > m_systemSize.x() * 0.5) pos[0] = pos[0] - m_systemSize.x();
    else if (pos.x() <= -m_systemSize.x() * 0.5) pos[0] = pos[0] + m_systemSize.x();
    if (pos.y() > m_systemSize.y() * 0.5) pos[1] = pos[1] - m_systemSize.y();
    else if (pos.y() <= -m_systemSize.y() * 0.5) pos[1] = pos[1] + m_systemSize.y();
    if (pos.z() > m_systemSize.z() * 0.5) pos[2] = pos[2] - m_systemSize.z();
    else if (pos.z() <= -m_systemSize.z() * 0.5) pos[2] = pos[2] + m_systemSize.z();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    // There is a small displacement of the system due to atoms starting in 0.
    // If not there is an ugly glitch where many atoms gets moved in the first step.
    for (int n = 0; n < (int) m_atoms.size(); n++) {
        if (m_atoms[n]->position.x() < - 0.2) m_atoms[n]->position[0] += m_systemSize.x();
        else if (m_atoms[n]->position.x() >= m_systemSize.x() - 0.2) m_atoms[n]->position[0] -= m_systemSize.x();
        if (m_atoms[n]->position.y() < -0.2) m_atoms[n]->position[1] += m_systemSize.y();
        else if (m_atoms[n]->position.y() >= m_systemSize.y() - 0.2) m_atoms[n]->position[1] -= m_systemSize.y();
        if (m_atoms[n]->position.z() < - 0.2) m_atoms[n]->position[2] += m_systemSize.z();
        else if (m_atoms[n]->position.z() >= m_systemSize.z() - 0.2) m_atoms[n]->position[2] -= m_systemSize.z();
    }
    // Make sure the right atoms are put in correct cells.
    m_celllist->emptyCells();
    m_celllist->calculateCellAtoms();
}

/*
 * Check this method, the momentum is zero to begin with.
 */
void System::removeMomentum() {
    // Initially, when the atoms are given random velocities, there is a non-zero net momentum. We don't want any drift in the system, so we need to remove it.
    
    // First we calculate the net momentum of the system.
    // We then subtract this sum from the velocity of each atom.
    
    velocityOfCM = vec3(0.0, 0.0, 0.0);
    velocityOfCMAfter = vec3(0.0, 0.0, 0.0);
    vec3 temp;
    double totalMass = 0;
    for (int n = 0; n < (int) m_atoms.size(); n++) {
        temp = m_atoms[n]->velocity * (m_atoms[n]->mass());
        velocityOfCM.add(temp);
        totalMass += m_atoms[n]->mass();
    }

    velocityOfCM = velocityOfCM / totalMass;
    for (int n = 0; n < (int) m_atoms.size(); n++) {
        m_atoms[n]->velocity = m_atoms[n]->velocity - velocityOfCM;
    }

    // Comparing the the velocity of the center of mass before and after removal.
    for (int n = 0; n < (int) m_atoms.size(); n++) {
        temp = m_atoms[n]->velocity * (m_atoms[n]->mass());
        velocityOfCMAfter.add(temp);
    }

    velocityOfCMAfter = velocityOfCMAfter / totalMass;
}

void System::resetForcesOnAllAtoms() {
    for (int i = 0; i < (int) m_celllist->listOfCells().size(); i++) {
        for (int j = 0; j < (int) m_celllist->listOfCells()[i]->atomsClose().size(); j++) {
            m_celllist->listOfCells()[i]->atomsClose()[j]->resetForce();
        }
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double iT, double rcut) {
    int totalNumberOfUnitCells = numberOfUnitCellsEachDimension * numberOfUnitCellsEachDimension * numberOfUnitCellsEachDimension;
    m_latticeConstant = latticeConstant;
    m_systemSize = vec3(latticeConstant * numberOfUnitCellsEachDimension, latticeConstant * numberOfUnitCellsEachDimension, latticeConstant * numberOfUnitCellsEachDimension);
    m_rcut = rcut;
    vec3 r1 = vec3(0.0, 0.0, 0.0);
    vec3 r2 = vec3(latticeConstant/2.0, latticeConstant/2.0, 0.0);
    vec3 r3 = vec3(0, latticeConstant/2.0, latticeConstant/2.0);
    vec3 r4 = vec3(latticeConstant/2.0, 0.0, latticeConstant/2.0);
    vector<vec3> R;
    vec3 temp;
    m_celllist = new CellList();
    m_celllist->setSystem(this);
    
    double initialTemperature = iT; // Dimensionless temperature.

    for (int i = 0; i < numberOfUnitCellsEachDimension; i++) {
        for (int j = 0; j < numberOfUnitCellsEachDimension; j++) {
            for (int k = 0; k < numberOfUnitCellsEachDimension; k++) {
                temp = vec3(i * latticeConstant, j * latticeConstant, k * latticeConstant);
                R.push_back(temp);
            }
        }
    }

    numberOfCellsX = int(m_systemSize.x() / m_rcut);
    numberOfCellsY = int(m_systemSize.y() / m_rcut);
    numberOfCellsZ = int(m_systemSize.z() / m_rcut);
    m_celllist->setNumberOfCells(numberOfCellsX, numberOfCellsY, numberOfCellsZ);
    for (int i = 0; i < numberOfCellsX; i++) {
        for (int j = 0; j < numberOfCellsY; j ++) {
            for (int k = 0; k < numberOfCellsZ; k++) {
                m_celllist->createCell();
            }
        }
    }

    // Creating atoms in unit cell formations.
    for (int n = 0; n < totalNumberOfUnitCells; n++) {
        Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom1->resetVelocityMaxwellian(initialTemperature);
        atom1->position = r1 + R[n];
        m_atoms.push_back(atom1);

        Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom2->resetVelocityMaxwellian(initialTemperature);
        atom2->position = r2 + R[n];
        m_atoms.push_back(atom2);

        Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom3->resetVelocityMaxwellian(initialTemperature);
        atom3->position = r3 + R[n];
        m_atoms.push_back(atom3);

        Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom4->resetVelocityMaxwellian(initialTemperature);
        atom4->position = r4 + R[n];
        m_atoms.push_back(atom4);
    }
    m_celllist->calculateCellAtoms();
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    if (m_oldForce) m_potential->calculateForcesOld(this);
    else m_potential->calculateForces(this);
}

void System::step(double dt)
{
    m_integrator->integrate(this, dt, m_thermostatOn);
    m_steps++;
    m_currentTime += dt;
}


void System::save(string filename)
{
    ofstream file(filename, ios::out | ios::binary);
    if (!file.is_open()) {
        cerr << "Could not open file " << filename << ". Aborting!" << endl;
        exit(1);
    }
    int numberOfPhaseSpaceCoordinates = 6 * m_atoms.size();
    double *phaseSpace = new double[numberOfPhaseSpaceCoordinates];

    int phaseSpaceCounter = 0;
    for (int i = 0; i < (int) m_atoms.size(); i++) {
        phaseSpace[phaseSpaceCounter++] = m_atoms[i]->position.x();
        phaseSpace[phaseSpaceCounter++] = m_atoms[i]->position.y();
        phaseSpace[phaseSpaceCounter++] = m_atoms[i]->position.z();
        phaseSpace[phaseSpaceCounter++] = m_atoms[i]->velocity.x();
        phaseSpace[phaseSpaceCounter++] = m_atoms[i]->velocity.y();
        phaseSpace[phaseSpaceCounter++] = m_atoms[i]->velocity.z();
    }

    int numberOfAtoms = m_atoms.size();
    file.write(reinterpret_cast<char*>(&numberOfAtoms), sizeof(int));
    file.write(reinterpret_cast<char*>(phaseSpace),
               numberOfPhaseSpaceCoordinates * sizeof(double));
    file.write(reinterpret_cast<char*>(&m_systemSize[0]), 3 * sizeof(double));

    file.close();
    delete phaseSpace;
}

void System::load(string filename)
{
    ifstream file(filename, ios::in | ios::binary);
    if (!file.is_open()) {
        cerr << "Could not open file " << filename << ". Aborting!" << endl;
        exit(1);
    }

    int numberOfAtoms;
    file.read(reinterpret_cast<char*>(&numberOfAtoms), sizeof(int));

    double *phaseSpace = new double[6 * numberOfAtoms];
    file.read(reinterpret_cast<char*>(phaseSpace), 6 * numberOfAtoms * sizeof(double));

    vec3 systemSize;
    file.read(reinterpret_cast<char*>(&systemSize[0]), 3 * sizeof(double));
    m_systemSize = systemSize;

    int phaseSpaceCounter = 0;
    double mass = 39.948;
    m_atoms.clear();
    for (int i = 0; i < numberOfAtoms; i++) {
        Atom *atom = new Atom(mass);
        atom->position = vec3(phaseSpace[phaseSpaceCounter++],
                              phaseSpace[phaseSpaceCounter++],
                              phaseSpace[phaseSpaceCounter++]);
        atom->velocity = vec3(phaseSpace[phaseSpaceCounter++],
                              phaseSpace[phaseSpaceCounter++],
                              phaseSpace[phaseSpaceCounter++]);
        m_atoms.push_back(atom);
    }

    file.close();
    delete phaseSpace;
}
