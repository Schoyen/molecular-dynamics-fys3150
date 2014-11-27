#include "lennardjones.h"
#include "../unitconverter.h"
#include <cmath>
#include <fstream>
#include "../celllist.h"
#include "../cell.h"

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    m_pressure = 0;
    double distanceBetweenAtoms = 0;
    double divisionOfSigmaAndDistance = 0;
    vec3 tempForce;
    vec3 distance;
    double expressionOfForce;
    int cx, cy, cz;

    //TODO: "Add the statistical property calculations to the cell lists";
    CellList *celllist = system->celllist();
    Cell *cell1;
    Cell *cell2;

    // And here we declare the winner of Big O.
    for (int i = 0; i < system->numberOfCellsX; i++) {
        //std::cout << "balle" << i << std::endl;
        for (int j = 0; j < system->numberOfCellsY; j++) {
            for (int k = 0; k < system->numberOfCellsZ; k++) {
                cell1 = celllist->getCell(i, j, k); // Current cell.

                for (int dx = i - 1; dx <= i + 1; dx++) {

                    cx = dx;
                    if (cx == system->numberOfCellsX) cx = 0;
                    else if (cx == -1) cx = system->numberOfCellsX - 1;

                    for (int dy = j - 1; dy <= j + 1; dy++) {

                        cy = dy;
                        if (cy == system->numberOfCellsY) cy = 0;
                        else if (cy == -1) cy = system->numberOfCellsY - 1;

                        for (int dz = k - 1; dz <= k + 1; dz++) {

                            cz = dz;
                            if (cz == system->numberOfCellsZ) cz = 0;
                            else if (cz == -1) cz = system->numberOfCellsZ - 1;

                            cell2 = celllist->getCell(cx, cy, cz);

                            if (cell1 == cell2) {
                                //std::cout << "Calculating locally" << std::endl;
                                //std::cout << i << "\t" << j << "\t" << k << "\n" << std::endl;
                                // Calulating force locally in a cell.

                                for (int m = 0; m < (int) cell1->atomsClose().size(); m++) {
                                    for (int n = m + 1; n < (int) cell1->atomsClose().size(); n++) {
                                        distance = cell1->atomsClose()[m]->position - cell1->atomsClose()[n]->position;
                                        distanceBetweenAtoms = distance.length();
                                        divisionOfSigmaAndDistance = m_sigma / distanceBetweenAtoms;

                                        m_potentialEnergy += 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                                        divisionOfSigmaAndDistance = m_sigma / system->rcut();
                                        m_potentialEnergy -= 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));

                                        expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12) / pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6) / pow(distanceBetweenAtoms, 8));
                                        tempForce = distance * expressionOfForce;
                                        cell1->atomsClose()[m]->force.add(tempForce);

                                        m_pressure += tempForce.dot(distance);

                                        tempForce = tempForce * (-1); // N3L
                                        cell1->atomsClose()[n]->force.add(tempForce);
                                    }
                                    // Calculate this in statisticssampler.
                                }
                            }
                            else {
                                // Calculate between cells.
                                //cell2 = celllist->getCell(cx, cy, cz); // Neighbouring cell.
                                
                                for (int m = 0; m < (int) cell1->atomsClose().size(); m++) {
                                    for (int n = 0; n < (int) cell2->atomsClose().size(); n++) {
                                        if (cell1->atomsClose()[m]->index >= cell2->atomsClose()[n]->index) {continue;}
                                        else {
                                            distance = cell1->atomsClose()[m]->position - cell2->atomsClose()[k]->position;
                                            // Write the minimumImageCriterion as a reference.
                                            system->minimumImageCriterion(distance);

                                            distanceBetweenAtoms = distance.length();
                                            divisionOfSigmaAndDistance = m_sigma / distanceBetweenAtoms;
                                            m_potentialEnergy += 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                                            expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12) / pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6) / pow(distanceBetweenAtoms, 8));
                                            tempForce = distance * expressionOfForce;
                                            cell1->atomsClose()[m]->force.add(tempForce);

                                            m_pressure += tempForce.dot(distance);

                                            tempForce = tempForce * (-1); // N3L
                                            cell2->atomsClose()[n]->force.add(tempForce);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}





void LennardJones::calculateForcesOld(System *system)
{
    m_potentialEnergy = 0;
    m_pressure = 0;
    double distanceBetweenAtoms = 0;
    double divisionOfSigmaAndDistance = 0;
    vec3 tempForce;
    vec3 distance;
    double expressionOfForce;

    for (int i = 0; i < (int) system->atoms().size(); i++) {
        for (int j = i + 1; j < (int) system->atoms().size(); j++) {
            distance = system->atoms()[i]->position - system->atoms()[j]->position;
            system->minimumImageCriterion(distance);
            if (distance.lengthSquared() > system->rcut() * system->rcut()) {
                continue;
            } else {
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
        }
    }
}
