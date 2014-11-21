#include "cell.h"
#include <cmath>

Cell::Cell()
{

}

Cell::~Cell() 
{

}

void Cell::clearList() {
    m_atomsClose.clear();
}

void Cell::addAtom(Atom *atom) {
    m_atomsClose.push_back(atom);
}

/*
 * This method is not "optimal".
 * Figure out a way to check whether or not an atom is inside a cell.
 */
bool Cell::isInCell(vec3 pos, double rcut)
{
    // This one needs to change.
}

/*
 * Method used to compute the forces between all atoms inside a cell locally.
 * Perhaps this method should be inline? It has been added here for readability reasons.
 */
double Cell::calculateLocally(double sigma, double epsilon, double rcut)
{
    vec3 distance;
    vec3 temp;
    double distanceBetweenAtoms;
    double divisionOfSigmaAndDistance;
    double potentialEnergy = 0;
    double expressionOfForce;
    for (int i = 0; i < (int) m_atomsClose.size(); i++) {
        for (int j = i + 1; j < (int) m_atomsClose.size(); j++) {
            distance = m_atomsClose[i]->position - m_atomsClose[j]->position;
            distanceBetweenAtoms = distance.length();
            divisionOfSigmaAndDistance = sigma / distanceBetweenAtoms;

            potentialEnergy += 4 * epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
            divisionOfSigmaAndDistance = sigma / rcut;
            potentialEnergy -= 4 * epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));

            expressionOfForce = 4 * epsilon * (12 * (pow(sigma, 12) / pow(distanceBetweenAtoms, 14)) - 6 * (pow(sigma, 6) / pow(distanceBetweenAtoms, 8)));
            temp = distance * expressionOfForce;
            m_atomsClose[i]->force.add(temp);
            temp = temp * (-1); // N3L
            m_atomsClose[j]->force.add(temp);
        }
    }
    calculatedLocally = true;
    return potentialEnergy;
}
