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
    vec3 temp;

    // TODO: Remember to update the celllists.
    // TODO: You must check whether or not the atoms in neighbouring cells are within r_cut.
    CellList *celllist = system->celllist();
    for (int i = 0; i < (int) celllist->listOfCells().size(); i++) {
        for (int j = i + 1; j < (int) celllist->listOfCells().size(); j++) {
            temp = celllist->listOfCells()[i]->position - celllist->listOfCells()[j]->position;
            //std::cout << tempTemp.lengthSquared() << std::endl;
            //std::cout << (pow(celllist->getrcut(), 2)) << std::endl;
            if (temp.lengthSquared() <= pow(celllist->getrcut(), 2)) {
                for (int k = 0; k < (int) celllist->listOfCells()[i]->atomsClose().size(); k++) {
                    // Calculation of force inside each cell.
                    for (int m = k + 1; m < (int) celllist->listOfCells()[i]->atomsClose().size(); m++) {
                        distance = celllist->listOfCells()[i]->atomsClose()[k]->position - celllist->listOfCells()[i]->atomsClose()[m]->position;

                        // Needed for minimum image criterion.
                        //x = distance.x();
                        //y = distance.y();
                        //z = distance.z();

                        // Do I need the minimum image criterion for the atoms locally?
                        /*
                        // Minimum image criterion.
                        if (x > xlen * 0.5) x = x - xlen;
                        else if (x < -xlen * 0.5) x = x + xlen;
                        if (y > ylen * 0.5) y = y - ylen;
                        else if (y < -ylen * 0.5) y = y + ylen;
                        if (z > zlen * 0.5) z = z - zlen;
                        else if (z < -zlen * 0.5) z = z + zlen;
                        */
                        
                        // Needed for minimum image criterion.
                        //distance = vec3(x, y, z);
                        
                        distanceBetweenAtoms = distance.length();
                        divisionOfSigmaAndDistance = m_sigma/distanceBetweenAtoms;
                        // Calculate the new potential energy.
                        // This one needs to be implemented properly.
                        m_potentialEnergy += 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                        divisionOfSigmaAndDistance = m_sigma/celllist->getrcut();
                        m_potentialEnergy -= 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                        expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12)/pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6)/pow(distanceBetweenAtoms, 8));
                        tempForce = distance * expressionOfForce;
                        celllist->listOfCells()[i]->atomsClose()[k]->force.add(tempForce);
                        tempForce = tempForce * (-1);
                        celllist->listOfCells()[i]->atomsClose()[m]->force.add(tempForce);
                    }

                    // Check whether of not he cells are within r_cut, at least for the potential.
                    // Calculation of force between each atom in neighbour cell[i] and cell[j].
                    for (int m = 0; m < (int) celllist->listOfCells()[j]->atomsClose().size(); m++) {
                        distance = celllist->listOfCells()[i]->atomsClose()[k]->position - celllist->listOfCells()[j]->atomsClose()[m]->position;
                        if (distance.lengthSquared() <= pow(celllist->getrcut(), 2)) {
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

                            m_potentialEnergy += 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                            divisionOfSigmaAndDistance = m_sigma/celllist->getrcut();
                            m_potentialEnergy -= 4 * m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
                            expressionOfForce = 4 * m_epsilon * (12 * pow(m_sigma, 12)/pow(distanceBetweenAtoms, 14) - 6 * pow(m_sigma, 6)/pow(distanceBetweenAtoms, 8));
                            tempForce = distance * expressionOfForce;
                            celllist->listOfCells()[i]->atomsClose()[k]->force.add(tempForce);
                            tempForce = tempForce * (-1);
                            celllist->listOfCells()[j]->atomsClose()[m]->force.add(tempForce);
                        }
                    }
                    m_kineticEnergy += 0.5 * celllist->listOfCells()[i]->atomsClose()[k]->mass() * celllist->listOfCells()[i]->atomsClose()[k]->velocity.lengthSquared();
                }
            }
        }
    }
    m_temperature = (2.0/3.0) * (m_kineticEnergy/((double) system->atoms().size() * 1));

    // TODO: Remember to compute temperature. The temperature is wrong.
    system->celllist()->emptyCells();
    system->celllist()->calculateCellAtoms();

    /*
    // Old force calculation.
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
