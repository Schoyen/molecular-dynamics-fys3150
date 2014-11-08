#include "statisticssampler.h"
#include <fstream>

StatisticsSampler::StatisticsSampler()
{
    sampledNetMomentum = false;
}

StatisticsSampler::~StatisticsSampler()
{

}

/*
 * The system seems to be working, but the temperature is 0 K and the
 * energies are fluctuating.
 */
void StatisticsSampler::sample(System *system, std::string filename)
{
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    if (!sampledNetMomentum) {
        sampleNetMomentum(system);
        sampledNetMomentum = true;
    }
    sampleTemperature(system);

    double temperatureInSI = UnitConverter::temperatureToSI(m_temperature);
    double kineticEnergyInSI = UnitConverter::energyToSI(m_kineticEnergy);
    double potentialEnergyInSI = UnitConverter::energyToSI(m_potentialEnergy);
    double netMomentumInSI = UnitConverter::velocityToSI(m_netMomentum);
    double netMomentumAfterInSI = UnitConverter::velocityToSI(m_netMomentumAfter);

    std::ofstream file(filename);
    if (file.is_open()) {
        file << "Initial net momentum of the system: " << netMomentumInSI << " m/s\n";
        file << "Net momentum of the system after removal: " << netMomentumAfterInSI << " m/s\n";
        file << "U = " << potentialEnergyInSI << " J\tK = " << kineticEnergyInSI << " J\tT = " << temperatureInSI << " K\n";
    } else std::cout << "Unable to write to file." << std::endl;
    file.close();
}

void StatisticsSampler::sampleKineticEnergy(System *system) {
    m_kineticEnergy = system->potential()->kineticEnergy();
}

void StatisticsSampler::samplePotentialEnergy(System *system) {
    m_potentialEnergy = system->potential()->potentialEnergy();
}

void StatisticsSampler::sampleNetMomentum(System *system) {
    // Value of the net momentum of the system before we remove it.
    m_netMomentum = system->getVelocityOfCM().length();
    m_netMomentumAfter = system->getVelocityOfCMAfter().length();
}

void StatisticsSampler::sampleTemperature(System *system) {
    m_temperature = system->potential()->temperature();
}
