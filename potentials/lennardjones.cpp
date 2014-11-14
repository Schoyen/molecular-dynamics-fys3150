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

/*
 * This does not look correct...
 */
void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    m_kineticEnergy = 0;
    double distanceBetweenAtoms = 0;
    double divisionOfSigmaAndDistance = 0;
    vec3 tempForce = vec3(0.0, 0.0, 0.0);
    vec3 distance;
    double expressionOfForce;
    double x, y, z;
    double xlen = system->systemSize().x();
    double ylen = system->systemSize().y();
    double zlen = system->systemSize().z();
    double counter = 0;
    vec3 temp;
    CellList *celllist = system->celllist();

    for (int i = 0; i < (int) celllist->listOfCells().size(); i++) {
        // Calculation of force inside each cell.
        // The local force seems to be computed properly.
        if (!celllist->listOfCells()[i]->calculatedLocally) {
            m_potentialEnergy += celllist->listOfCells()[i]->calculateLocally(m_sigma, m_epsilon, celllist->getrcut());
        }


        for (int j = i + 1; j < (int) celllist->listOfCells().size(); j++) {

            temp = celllist->listOfCells()[i]->position - celllist->listOfCells()[j]->position;
            x = temp.x();
            y = temp.y();
            z = temp.z();

            // Minimum image criterion for the cells.
            // This is needed for calculation of cells the are along the edges.
            if (x > xlen * 0.5) x = x - xlen;
            else if (x < -xlen * 0.5) x = x + xlen;
            if (y > ylen * 0.5) y = y - ylen;
            else if (y < -ylen * 0.5) y = y + ylen;
            if (z > zlen * 0.5) z = z - zlen;
            else if (z < -zlen * 0.5) z = z + zlen;
            temp = vec3(x, y, z);

            // This isn't working.
            // Find a proper test.
            if (temp.x() <= celllist->listOfCells()[j]->position.x() + celllist->getrcut() &&
                temp.y() <= celllist->listOfCells()[j]->position.y() + celllist->getrcut() &&
                temp.z() <= celllist->listOfCells()[j]->position.z() + celllist->getrcut()) {
                counter++;
                for (int k = 0; k < (int) celllist->listOfCells()[i]->atomsClose().size(); k++) {
                    //std::cout << "Calculating force between atoms" << std::endl;
                    // Check whether of not the cells are within r_cut, at least for the potential.
                    // Calculation of force between each atom in neighbour cell[i] and cell[j].
                    for (int m = 0; m < (int) celllist->listOfCells()[j]->atomsClose().size(); m++) {
                        distance = celllist->listOfCells()[i]->atomsClose()[k]->position - celllist->listOfCells()[j]->atomsClose()[m]->position;
                        x = distance.x();
                        y = distance.y();
                        z = distance.z();

                        // Minimum image criterion.
                        if (x > xlen * 0.5) x = x - xlen;
                        else if (x < -xlen * 0.5) x = x + xlen;
                        if (y > ylen * 0.5) y = y - ylen;
                        else if (y < -ylen * 0.5) y = y + ylen;
                        if (z > zlen * 0.5) z = z - zlen;
                        else if (z < -zlen * 0.5) z = z + zlen;
                        
                        distance = vec3(x, y, z);
                        distanceBetweenAtoms = distance.length();
                        divisionOfSigmaAndDistance = m_sigma/distanceBetweenAtoms;

                        // Calculating the potential according to the formula.
                        if (distance.lengthSquared() <= 3 * pow(celllist->getrcut(), 2)) {
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
        std::cout << counter << std::endl;
        counter = 0;
    }
    m_temperature = (2.0/3.0) * (m_kineticEnergy/((double) system->atoms().size() * 1));








    // Needed for timing of the methods.
    // Old force calculation.
    /*
    for (int i = 0; i < (int) system->atoms().size(); i++) {
        for (int j = i + 1; j < (int) system->atoms().size(); j++) {
            distance = system->atoms()[i]->position - system->atoms()[j]->position;
            x = distance.x();
            y = distance.y();
            z = distance.z();

            // Minimum image criterion.
            if (x > xlen * 0.5) x = x - xlen;
            else if (x < -xlen * 0.5) x = x + xlen;
            if (y > ylen * 0.5) y = y - ylen;
            else if (y < -ylen * 0.5) y = y + ylen;
            if (z > zlen * 0.5) z = z - zlen;
            else if (z < -zlen * 0.5) z = z + zlen;
            distance = vec3(x, y, z);
            distanceBetweenAtoms = distance.length();
            divisionOfSigmaAndDistance = m_sigma/distanceBetweenAtoms;
            m_potentialEnergy += 4 * m_epsilon * (pow(distanceBetweenAtoms, 12) - pow(divisionOfSigmaAndDistance, 6));
            expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12)/pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6)/pow(distanceBetweenAtoms, 8));
            tempForce = distance*expressionOfForce;
            system->atoms()[i]->force.add(tempForce);
            tempForce = tempForce * (-1);
            system->atoms()[j]->force.add(tempForce);
        }
        m_kineticEnergy += 0.5 * system->atoms()[i]->mass() * system->atoms()[i]->velocity.lengthSquared();
    }
    m_temperature = (2.0/3.0) * (m_kineticEnergy/((double) system->atoms().size() * 1));
    */
}
