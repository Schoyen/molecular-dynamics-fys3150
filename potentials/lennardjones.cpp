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
    int counter;
    int cx, cy, cz;

    //TODO: "Add the statistical property calculations to the cell lists";
    CellList *celllist = system->celllist();
    Cell *cell1;
    Cell *cell2;

    // And here we declare the winner of Big O.
    for (int i = 0; i < system->numberOfCellsX; i++) {
        for (int j = 0; j < system->numberOfCellsY; j++) {
            for (int k = 0; k < system->numberOfCellsZ; k++) {
                cell1 = celllist->getCell(i, j, k); // Current cell.
                counter = 0;

                // An attempt to improve speed.
                // We are only calculating neigbour cells in front of us.
                // Due to the minimum image criterion this should work.
                for (int dx = i; dx <= i + 1; dx++) {

                    cx = dx;
                    if (cx == system->numberOfCellsX) cx = 0;
                    // This one should only be necessary if we calculate for all cells.
                    //else if (cx == -1) cx = system->numberOfCellsX - 1;

                    for (int dy = j; dy <= j + 1; dy++) {

                        cy = dy;
                        if (cy == system->numberOfCellsY) cy = 0;
                        //else if (cy == -1) cy = system->numberOfCellsY - 1;

                        for (int dz = k; dz <= k + 1; dz++) {

                            cz = dz;
                            if (cz == system->numberOfCellsZ) cz = 0;
                            //else if (cz == -1) cz = system->numberOfCellsZ - 1;

                            if (counter == 0) {
                                std::cout << i << "\t" << j << "\t" << k << "\n" << std::endl;
                                // Calulating force locally in a cell.

                                for (int k = 0; i < (int) cell1->atomsClose().size(); i++) {
                                    for (int m = k + 1; m < (int) cell1->atomsClose().size(); m++) {
                                        distance = cell1->atomsClose()[k]->position - cell1->atomsClose()[m]->position;
                                        distanceBetweenAtoms = distance.length();
                                        divisionOfSigmaAndDistance = m_sigma / distanceBetweenAtoms;

                                        m_potentialEnergy += 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                                        divisionOfSigmaAndDistance = m_sigma / celllist->getrcut();
                                        m_potentialEnergy -= 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));

                                        expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12) / pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6) / pow(distanceBetweenAtoms, 8));
                                        tempForce = distance * expressionOfForce;
                                        cell1->atomsClose()[k]->force.add(tempForce);

                                        m_pressure += tempForce.dot(distance);

                                        tempForce = tempForce * (-1); // N3L
                                        cell1->atomsClose()[m]->force.add(tempForce);
                                    }
                                    m_kineticEnergy += 0.5 * cell1->atomsClose()[k]->mass() * cell1->atomsClose()[k]->velocity.lengthSquared();
                                }
                            }
                            else {
                                // Calculate between cells.
                                cell2 = celllist->getCell(cx, cy, cz); // Neighbouring cell.
                                std::cout << cx << "\t" << cy << "\t" << cz << std::endl;
                                
                                for (int m = 0; m < (int) cell1->atomsClose().size(); m++) {
                                    for (int k = 0; k < (int) cell2->atomsClose().size(); k++) {
                                        distance = cell1->atomsClose()[m]->position - cell2->atomsClose()[k]->position;
                                        distance = system->minimumImageCriterion(distance);

                                        distanceBetweenAtoms = distance.length();
                                        divisionOfSigmaAndDistance = m_sigma / distanceBetweenAtoms;
                                        m_potentialEnergy += 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                                        expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12) / pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6) / pow(distanceBetweenAtoms, 8));
                                        tempForce = distance * expressionOfForce;
                                        cell1->atomsClose()[m]->force.add(tempForce);

                                        m_pressure += tempForce.dot(distance);

                                        tempForce = tempForce * (-1); // N3L
                                        cell2->atomsClose()[k]->force.add(tempForce);
                                    }
                                    m_kineticEnergy += 0.5 * cell1->atomsClose()[m]->mass() * cell1->atomsClose()[m]->velocity.lengthSquared();
                                }
                            }
                            counter++;
                        }
                    }
                }
            }
        }
    }
    m_temperature = (2.0/3.0) * (m_kineticEnergy / ((double) system->atoms().size() * 1));
    m_pressure = 1.0 / (3.0 * system->systemSize().x() * system->systemSize().y() * system->systemSize().z() * (m_pressure / ((double) system->atoms().size()))) + m_numberDensity * 1 * m_temperature;
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

    for (int i = 0; i < (int) system->atoms().size(); i++) {
        for (int j = i + 1; j < (int) system->atoms().size(); j++) {
            distance = system->atoms()[i]->position - system->atoms()[j]->position;
            distance = system->minimumImageCriterion(distance);

            distanceBetweenAtoms = distance.length();
            divisionOfSigmaAndDistance = m_sigma/distanceBetweenAtoms;
            m_potentialEnergy += 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
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
