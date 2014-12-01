#include <iostream>
#include <cstdlib>
#include <chrono>
#include "io.h"
#include "system.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "berendsen.h"
#include "potentials/lennardjones.h"
#include "integrators/velocityverlet.h"

using namespace std;
using namespace chrono;
/*
 * List of things that need to be sent to the main-method.
 * double dt
 * int numberOfFCCLattices
 * int cellSize
 * double initialTemperature
 * double tbath
 * double relaxationTime
 * int timestep
 * int timestepStartThermostat
 * int timestepEndThermostat
 * double latticeConstant
 * double rcut
 * bool oldForceCalculation
 *
 * The program should save the state of the system after it has been completed and automatically check if there is an existing loaded state.
 */
int main(int argc, char *argv[])
{
    if (argc < 12) {
        cout << "\n==================================" << endl;
        cout << "Not enough command line arguments." << endl;
        cout << "Usage: " << argv[0] << " 1 2 3 4 5 6 7 8 9 10 11 12 13\n"
             << "1: double dt\n"
             << "2: int number of FCC lattices\n"
             << "3: double initial temperature\n"
             << "4: double t_bath\n"
             << "5: double relaxation time\n"
             << "6: int number of timesteps\n"
             << "7: int number of timesteps before turning on the thermostat\n"
             << "8: int number of timesteps before turning off the thermostat\n"
             << "9: double latticConstant\n"
             << "10: double rcut\n"
             << "11: bool 1='true' for old force calculation and 0 (or anything) ='false' for cell lists"
             << endl;
        cout << "==================================\n" << endl;
        exit(1);
    }

    double dt = UnitConverter::timeFromSI(atof(argv[1]));
    int numberOfFCCLattices = atoi(argv[2]);
    double initialTemperature = atof(argv[3]);
    double tbath = atof(argv[4]);
    double relaxationTime = atof(argv[5]);
    int timesteps = atoi(argv[6]);
    int timestepStartThermostat = atoi(argv[7]);
    int timestepEndThermostat = atoi(argv[8]);
    double latticeConstant = atof(argv[9]);
    double rcut = atof(argv[10]);
    bool oldForceCalculation;
    if (atoi(argv[11]) == 1) oldForceCalculation = true;
    else oldForceCalculation = false;
    double sigma = 3.405; // From assignment text.
    double epsilon = 1.0;

    IO *movie = new IO();
    movie->open("movie.xyz");
    StatisticsSampler *statisticsSampler = new StatisticsSampler();
    System system;
    system.createFCCLattice(numberOfFCCLattices, UnitConverter::lengthFromAngstroms(latticeConstant),
                            UnitConverter::temperatureFromSI(initialTemperature), rcut);
    system.setPotential(new LennardJones(sigma, epsilon));
    system.setIntegrator(new VelocityVerlet());
    system.setThermostat(new BerendsenThermostat(UnitConverter::temperatureFromSI(tbath),
                                                 relaxationTime));
    system.removeMomentum();
    system.setForceCalculation(oldForceCalculation);
    system.setThermostatOn(false);
    statisticsSampler->createFiles();

    // Starting clock.
    auto start = high_resolution_clock::now();
    for (int i = 0; i < timesteps; i++) {
        movie->saveState(&system);
        system.temperature = statisticsSampler->temperature();

        // timestepStartThermostat == -1 if the thermostat should be off for all calculations.
        // These are not working properly.
        //if (i == timestepStartThermostat) system.setThermostatOn(true);
        //if (i == timestepEndThermostat) system.setThermostatOn(false);
        system.step(dt);
        statisticsSampler->sampleKineticEnergySquared(&system);
        statisticsSampler->sampleTotalKineticEnergy(&system);

        // Sampling every 100'th step.
        if (i % 100 == 0) statisticsSampler->sample(&system, i);
    }
    statisticsSampler->sampleHeatCapacity(&system);

    // Ending clock.
    auto finish = high_resolution_clock::now();
    double time = duration_cast<seconds>(finish - start).count();
    double kiloAtomPerTimestep = 4 * numberOfFCCLattices * numberOfFCCLattices * numberOfFCCLattices * timesteps / time;
    string filename = "build/DATA/kiloAtomTimestep.txt";
    ofstream file;
    file.open(filename);
    if (!file.is_open()) {
        cerr << "Unable to write to build/DATA/kiloAtomTimestep.txt" << endl;
        exit(1);
    }
    file << kiloAtomPerTimestep << "\n";
    file.close();

    movie->saveState(&system);
    statisticsSampler->closeFiles();
    movie->close();

    return 0;
}
