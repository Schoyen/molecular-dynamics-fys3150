#include "statisticssampler.h"

// Calculate kinetic energy, pressure, temperature and density in here.

using namespace std;

StatisticsSampler::StatisticsSampler()
{
    m_sampledNetMomentum = false;
    m_sampledNumberDensity = false;
}

StatisticsSampler::~StatisticsSampler()
{

}


void StatisticsSampler::createFiles()
{
    m_kineticEnergyFile.open("build/DATA/kineticEnergy.txt");
    m_potentialEnergyFile.open("build/DATA/potentialEnergy.txt");
    m_temperatureFile.open("build/DATA/temperature.txt");
    m_netMomentumFile.open("build/DATA/netMomentum.txt");
    m_numberDensityFile.open("build/DATA/numberDensity.txt");
    m_pressureFile.open("build/DATA/pressure.txt");
    // Consider calculating this in a python program when plotting.
    // Check whether or not you are supposed to use energy divided by temperature or some other way.
    m_heatCapacityFile.open("build/DATA/heatCapacity.txt");
}


void StatisticsSampler::closeFiles()
{
    m_kineticEnergyFile.close();
    m_potentialEnergyFile.close();
    m_temperatureFile.close();
    m_netMomentumFile.close();
    m_numberDensityFile.close();
    m_pressureFile.close();
    m_heatCapacityFile.close();
}

/*
 * Figure out trouble with the net momentum. It is zero to begin with.
 */
void StatisticsSampler::sample(System *system, int timestep)
{
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    samplePressure(system);
    if (!m_sampledNetMomentum) {
        sampleNetMomentum(system);
        m_sampledNetMomentum = true;
        double netMomentumInSI = UnitConverter::velocityToSI(m_netMomentum);
        double netMomentumAfterInSI = UnitConverter::velocityToSI(m_netMomentumAfter);
        if (!m_netMomentumFile.is_open()) {
            cerr << "Unable to write to build/DATA/netMomentum.txt" << endl;
            exit(1);
        }
        m_netMomentumFile << "Net momentum before removal: " << to_string(netMomentumInSI) << "\n";
        m_netMomentumFile << "Net momentum after removal: " << to_string(netMomentumAfterInSI);
    }

    if (!m_sampledNumberDensity) {
        sampleNumberDensity(system);
        m_sampledNumberDensity = true;
        if (!m_numberDensityFile.is_open()) {
            cerr << "Unable to write to build/DATA/numberDensity.txt" << endl;
            exit(1);
        }

        m_numberDensityFile << m_numberDensity << "\n";
    }


    double temperatureInSI = UnitConverter::temperatureToSI(m_temperature);
    double kineticEnergyInSI = UnitConverter::energyToEv(m_kineticEnergy);
    double potentialEnergyInSI = UnitConverter::energyToEv(m_potentialEnergy);
    double pressureInSI = UnitConverter::pressureToSI(m_pressure);

    if (!m_kineticEnergyFile.is_open()) {
        cerr << "Unable to write to build/DATA/kineticEnergy.txt" << endl;
        exit(1);
    }
    
    m_kineticEnergyFile << kineticEnergyInSI << "\n";

    if (!m_potentialEnergyFile.is_open()) {
        cerr << "Unable to write to build/DATA/potentialEnergy.txt" << endl;
        exit(1);
    }
    
    m_potentialEnergyFile << potentialEnergyInSI << "\n";

    if (!m_temperatureFile.is_open()) {
        cerr << "Unable to write to build/DATA/temperature.txt" << endl;
        exit(1);
    }
    
    m_temperatureFile << temperatureInSI << "\n";


    if (!m_pressureFile.is_open()) {
        cerr << "Unable to write to build/DATA/pressure.txt" << endl;
        exit(1);
    }

    m_pressureFile << pressureInSI << "\n";

    // Add the extra conversions for the rest of the files.
}

void StatisticsSampler::sampleKineticEnergy(System *system)
{
    m_kineticEnergy = 0;
    for (int i = 0; i < (int) system->atoms().size(); i++) {
        m_kineticEnergy += 0.5 * system->atoms()[i]->mass() * system->atoms()[i]->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System *system)
{
    m_potentialEnergy = system->potential()->potentialEnergy();
}

void StatisticsSampler::sampleNetMomentum(System *system)
{
    // Value of the net momentum of the system before we remove it.
    // Both are zero. Make sure the function is working properly.
    m_netMomentum = system->getVelocityOfCM().length();
    m_netMomentumAfter = system->getVelocityOfCMAfter().length();
}

void StatisticsSampler::sampleTemperature(System *system)
{
    m_temperature = (2.0/3.0) * (m_kineticEnergy / (system->atoms().size() * 1.0));
}

void StatisticsSampler::sampleNumberDensity(System *system)
{
    m_numberDensity = system->atoms().size() / ((double) system->systemSize().x() * system->systemSize().y() * system->systemSize().z());
}

void StatisticsSampler::samplePressure(System *system)
{
    m_pressure = system->potential()->pressure();
    m_pressure = 1.0 / (3.0 * system->systemSize().x() * system->systemSize().y() * system->systemSize().z()) * (m_pressure / ((double) system->atoms().size())) + m_numberDensity * 1.0 * m_temperature;
}
