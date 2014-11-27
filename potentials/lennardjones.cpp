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

    // TODO:
    // Declare variables in loop.
    // Remove pow.
    // Add booleans in regular expressions. TEST.
    // Use add and multiply.
    double distanceBetweenAtoms = 0;
    double divisionOfSigmaAndDistance = 0;
    vec3 tempForce;
    double expressionOfForce;
    int cx, cy, cz;

    CellList *celllist = system->celllist();
    Cell *cell1;
    Cell *cell2;

    // And here we declare the winner of Big O.
    for (int i = 0; i < system->numberOfCellsX; i++) {
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
                            // Calculate between cells.
                                
                            for (int m = 0; m < (int) cell1->atomsClose().size(); m++) {
                                for (int n = 0; n < (int) cell2->atomsClose().size(); n++) {
                                    if (cell1->atomsClose()[m]->index >= cell2->atomsClose()[n]->index) {continue;}
                                    else {
                                        vec3 distance = cell1->atomsClose()[m]->position;
                                        distance.addAndMultiply(cell2->atomsClose()[n]->position, -1);
                                        system->minimumImageCriterion(distance);
                                        
                                        if (distance.lengthSquared() > system->rcut() * system->rcut()) continue;

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
