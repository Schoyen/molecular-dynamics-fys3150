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
    vec3 tempForce = 0; // N3L
    vec3 distance;
    for (int i = 0; i < system->atoms().size(); i ++) {
        for (int j = i + 1; j < system->atoms().size(); j++) {
            distance = system->atoms[i] - system->atoms[j];
            distanceBetweenAtoms = distance.length();
            divisionOfSigmaAndDistance = m_sigma/distanceBetweenAtoms;
            m_potentialEnergy += 4*m_epsilon(pow(divisionOfSigmaAndDistance, 12) - pow(divisionOfSigmaAndDistance, 6));
            // You were here!
            tempForce = 
        }
    }
}
