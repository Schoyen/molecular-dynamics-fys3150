#pragma once
#include <vector>
#include "atom.h"
#include "math/vec3.h"
#include "potentials/potential.h"
#include "berendsen.h"
#include "celllist.h"

class Integrator; class Potential; class CellList;
using std::vector;
using CompPhys::vec3;
using namespace std;

class System
{
private:
    vec3 m_systemSize;
    vector<Atom*> m_atoms;
    Potential *m_potential;
    Integrator *m_integrator;
    double m_currentTime;
    double m_rcut;
    double m_latticeConstant;
    int m_steps;
    vec3 velocityOfCM;
    vec3 velocityOfCMAfter;
    CellList *m_celllist;
    BerendsenThermostat *m_thermostat;
    bool m_oldForce;
    bool m_thermostatOn;

public:
    int numberOfCellsX;
    int numberOfCellsY;
    int numberOfCellsZ;
    double temperature;

    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double iT, double rcut);
    void applyPeriodicBoundaryConditions();
    void minimumImageCriterion(vec3 &pos);
    void removeMomentum();
    void calculateForces();
    void step(double dt);
    void save(string filename);
    void load(string filename);

    // Setters and getters
    vector<Atom *> &atoms() { return m_atoms; }
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    void setForceCalculation(bool oldForce) {m_oldForce = oldForce;}
    void setThermostatOn(bool thermostatOn) {m_thermostatOn = thermostatOn;}
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double currentTime() { return m_currentTime; }
    double rcut() {return m_rcut;}
    void setCurrentTime(double currentTime) { m_currentTime = currentTime; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    vec3 getVelocityOfCM() { return velocityOfCM; }
    vec3 getVelocityOfCMAfter() { return velocityOfCMAfter; }
    CellList *celllist() { return m_celllist; }
    void setThermostat(BerendsenThermostat *thermostat) {m_thermostat = thermostat;}
    BerendsenThermostat *berendsen() {return m_thermostat;}
    double latticeConstant() {return m_latticeConstant;}
};
