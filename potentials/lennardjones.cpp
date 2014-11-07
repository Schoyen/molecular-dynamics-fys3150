#include "lennardjones.h"
#include "../unitconverter.h"
#include <cmath>
#include <fstream>

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    m_kineticEnergy = 0;
    double distanceBetweenAtoms = 0;
    double divisionOfSigmaAndDistance = 0;
    vec3 tempForce = vec3(0.0, 0.0, 0.0); // N3L
    vec3 distance;
    double expressionOfForce;
    double x, y, z;
    double xlen = system->systemSize().x();
    double ylen = system->systemSize().y();
    double zlen = system->systemSize().z();
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
}
