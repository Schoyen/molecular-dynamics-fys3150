#include "lennardjones.h"
#include <cmath>

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    double distanceBetweenAtoms = 0;
    double divisionOfSigmaAndDistance = 0;
    double forceX;
    double forceY;
    double forceZ;
    vec3 tempForce = vec3(0.0, 0.0, 0.0); // N3L
    vec3 distance;
    for (int i = 0; i < (int) system->atoms().size(); i ++) {
        for (int j = i + 1; j < (int) system->atoms().size(); j++) {
            distance = system->atoms()[i]->position - system->atoms()[j]->position;
            distanceBetweenAtoms = distance.length();
            divisionOfSigmaAndDistance = m_sigma/distanceBetweenAtoms;
            m_potentialEnergy += 4*m_epsilon * (pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
            forceX = 4*m_epsilon*distance.x()*((12*pow(m_sigma, 12))/pow(distanceBetweenAtoms, 14) - (6*pow(m_sigma, 6))/pow(distanceBetweenAtoms, 8));
            forceY = 4*m_epsilon*distance.y()*((12*pow(m_sigma, 12))/pow(distanceBetweenAtoms, 14) - (6*pow(m_sigma, 6))/pow(distanceBetweenAtoms, 8));
            forceZ = 4*m_epsilon*distance.z()*((12*pow(m_sigma, 12))/pow(distanceBetweenAtoms, 14) - (6*pow(m_sigma, 6))/pow(distanceBetweenAtoms, 8));
            tempForce = vec3(forceX, forceY, forceZ);
            system->atoms()[i]->force.add(tempForce);
            tempForce = vec3(-forceX, -forceY, -forceZ);
            system->atoms()[j]->force.add(tempForce);
        }
    }
}
