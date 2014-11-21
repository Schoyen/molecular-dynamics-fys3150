#include "lennardjones.h"
#include "../unitconverter.h"
#include <cmath>
#include <fstream>
#include "../celllist.h"
#include "../cell.h"

LennardJones::LennardJones(double sigma, double epsilon, BerendsenThermostat *berendsen) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{
    m_berendsen = berendsen;
}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    m_kineticEnergy = 0;
    m_numberDensity = system->atoms().size() / (system->systemSize().x() * system->systemSize().y() * system->systemSize().z());
    double distanceBetweenAtoms = 0;
    double divisionOfSigmaAndDistance = 0;
    vec3 tempForce;
    vec3 distance;
    double expressionOfForce;
    vec3 temp;

    //TODO: "Add the statistical property calculations to the cell lists";
    CellList *celllist = system->celllist();
    Cell *cell1;
    Cell *cell2;

    for (int i = 0; i < system->numberOfCellsX; i++) {
        for (int j = 0; j < system->numberOfCellsY; j++) {
            for (int k = 0; k < system->numberOfCellsZ; k++) {
                cell1 = celllist->getCell(i, j, k);

                for (int dx = i - 1; dx <= i + 1; dx++) {
                    for (int dy = j - 1; dy <= j + 1; dy++) {
                        for (int dz = k - 1; dz <= k + 1; dz++) {
                            
                        }
                    }
                }
            }
        }
    }
}





void LennardJones::calculateForcesOld(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    m_kineticEnergy = 0;
    m_numberDensity = system->atoms().size() / (system->systemSize().x() * system->systemSize().y() * system->systemSize().z());
    double distanceBetweenAtoms = 0;
    double divisionOfSigmaAndDistance = 0;
    vec3 tempForce = vec3(0.0, 0.0, 0.0);
    vec3 distance;
    double expressionOfForce;
    vec3 temp;

    for (int i = 0; i < (int) system->atoms().size(); i++) {
        for (int j = i + 1; j < (int) system->atoms().size(); j++) {
            distance = system->atoms()[i]->position - system->atoms()[j]->position;
            distance = system->minimumImageCriterion(distance);

            distanceBetweenAtoms = distance.length();
            divisionOfSigmaAndDistance = m_sigma/distanceBetweenAtoms;
            m_potentialEnergy += 4 * m_epsilon * (pow(distanceBetweenAtoms, 12) - pow(divisionOfSigmaAndDistance, 6));
            expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12)/pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6)/pow(distanceBetweenAtoms, 8));
            tempForce = distance*expressionOfForce;
            system->atoms()[i]->force.add(tempForce);

            m_pressure += tempForce.dot(distance);

            tempForce = tempForce * (-1);
            system->atoms()[j]->force.add(tempForce);
        }
        m_kineticEnergy += 0.5 * system->atoms()[i]->mass() * system->atoms()[i]->velocity.lengthSquared();
    }
    m_temperature = (2.0/3.0) * (m_kineticEnergy/((double) system->atoms().size() * 1));

    m_pressure = 1.0 / (3.0 * system->systemSize().x() * system->systemSize().y() * system->systemSize().z()) * (m_pressure / ((double) system->atoms().size())) + m_numberDensity * 1 * m_temperature;
}
