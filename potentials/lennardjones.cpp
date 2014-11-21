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
    vec3 tempForce = vec3(0.0, 0.0, 0.0);
    vec3 distance;
    double expressionOfForce;
    double counter = 0;
    vec3 temp;

    //TODO: "Add the statistical property calculations to the cell lists";
    CellList *celllist = system->celllist();

    for (int i = 0; i < (int) celllist->listOfCells().size(); i++) {
        // Calculation of force inside each cell.
        // The local force seems to be computed properly.
        if (!celllist->listOfCells()[i]->calculatedLocally) {
            m_potentialEnergy += celllist->listOfCells()[i]->calculateLocally(m_sigma, m_epsilon, celllist->getrcut());
        }


        for (int j = i + 1; j < (int) celllist->listOfCells().size(); j++) {

            temp = celllist->listOfCells()[i]->position - celllist->listOfCells()[j]->position;
            temp = system->minimumImageCriterion(temp);

            // This isn't working.
            // Find a proper test.
            if (celllist->listOfCells()[j]->isInCell(temp, celllist->getrcut())) {
                counter++;
                for (int k = 0; k < (int) celllist->listOfCells()[i]->atomsClose().size(); k++) {
                    //std::cout << "Calculating force between atoms" << std::endl;
                    for (int m = 0; m < (int) celllist->listOfCells()[j]->atomsClose().size(); m++) {
                        distance = celllist->listOfCells()[i]->atomsClose()[k]->position - celllist->listOfCells()[j]->atomsClose()[m]->position;
                        distance = system->minimumImageCriterion(distance);
                        distanceBetweenAtoms = distance.length();
                        divisionOfSigmaAndDistance = m_sigma/distanceBetweenAtoms;

                        // Calculating the potential according to the formula.
                        if (celllist->listOfCells()[j]->isInCell(distance, celllist->getrcut())) {
                            m_potentialEnergy += 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                            divisionOfSigmaAndDistance = m_sigma/celllist->getrcut();
                            m_potentialEnergy -= 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                        }

                        expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12)/pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6)/pow(distanceBetweenAtoms, 8));
                        tempForce = distance * expressionOfForce;
                        celllist->listOfCells()[i]->atomsClose()[k]->force.add(tempForce);
                        tempForce = tempForce * (-1);
                        celllist->listOfCells()[j]->atomsClose()[m]->force.add(tempForce);
                    }
                    m_kineticEnergy += 0.5 * celllist->listOfCells()[i]->atomsClose()[k]->mass() * celllist->listOfCells()[i]->atomsClose()[k]->velocity.lengthSquared();
                }
            }
        }
        //std::cout << counter << std::endl;
        counter = 0;
    }
    m_temperature = (2.0/3.0) * (m_kineticEnergy/((double) system->atoms().size() * 1));








    // Needed for timing of the methods.
    // Old force calculation.
    // TODO: Write the old force calculation in a different method.
    /*
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
    */

    m_pressure = 1.0 / (3.0 * system->systemSize().x() * system->systemSize().y() * system->systemSize().z()) * (m_pressure / ((double) system->atoms().size())) + m_numberDensity * 1 * m_temperature;

}
