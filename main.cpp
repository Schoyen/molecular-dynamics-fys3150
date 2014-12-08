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
 * int timestepWithThermostat
 * double latticeConstant
 * double rcut
 * bool oldForceCalculation
 *
 * The program should save the state of the system after it has been completed and automatically check if there is an existing loaded state.
 */
int main(int argc, char *argv[])
{
    if (argc < 11) {
        cout << "\n==================================" << endl;
        cout << "Not enough command line arguments." << endl;
        cout << "Usage: " << argv[0] << " 1 2 3 4 5 6 7 8 9 10 \n"
             << "1: double dt\n"
             << "2: int number of FCC lattices\n"
             << "3: double initial temperature\n"
             << "4: double t_bath\n"
             << "5: double relaxation time\n"
             << "6: int number of timesteps\n"
             << "7: int number of timesteps with thermostat\n"
             << "8: double latticConstant\n"
             << "9: double rcut\n"
             << "10: bool 1='true' for old force calculation and 0 (or anything) ='false' for cell lists"
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
    int timestepWithThermostat = atoi(argv[7]);
    double latticeConstant = atof(argv[8]);
    double rcut = atof(argv[9]);
    bool oldForceCalculation;
    if (atoi(argv[10]) == 1) oldForceCalculation = true;
    else oldForceCalculation = false;
    double sigma = 3.405; // From assignment text.
    double epsilon = 1.0; // This should maybe be 119.8 K.

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
    string filename = "thermostat.txt";
    // Starting clock.
    if (argc == 11) {
        auto start = high_resolution_clock::now();
        for (int i = 0; i < timesteps; i++) {
            movie->saveState(&system);
            statisticsSampler->sample(&system, i);
            system.temperature = statisticsSampler->temperature();
            system.step(dt);
            statisticsSampler->sampleKineticEnergySquared(&system);
            statisticsSampler->sampleTotalKineticEnergy(&system);
            statisticsSampler->sampleKineticEnergy(&system);
            statisticsSampler->sampleTemperature(&system);
        }
        statisticsSampler->sampleHeatCapacity(&system);

        // Ending clock.
        auto finish = high_resolution_clock::now();
        double time = duration_cast<seconds>(finish - start).count();
        // Divided by a thousand?
        double kiloAtomPerTimestep = 4 * numberOfFCCLattices * numberOfFCCLattices * numberOfFCCLattices * timesteps / time * 0.001;
        int numberOfAtoms = (int) 4 * numberOfFCCLattices * numberOfFCCLattices * numberOfFCCLattices;
        string filename = "build/DATA/kiloAtomTimestep-" + to_string(numberOfAtoms) + "-" + to_string(timesteps) + "-" + argv[10] + ".txt";
        ofstream file;
        file.open(filename);
        if (!file.is_open()) {
            cerr << "Unable to write to " << filename << endl;
            exit(1);
        }
        file << kiloAtomPerTimestep << "\t" << time << "\t" << numberOfAtoms << "\n";
        file.close();
    } else if (atoi(argv[11]) == 1) {
        // Running the thermostat and saving the state after finishing.
        system.setThermostatOn(true);
        for (int i = 0; i < timestepWithThermostat; i++) {
            movie->saveState(&system);
            statisticsSampler->sample(&system, i);
            system.temperature = statisticsSampler->temperature();
            system.step(dt);
            statisticsSampler->sampleKineticEnergy(&system);
            statisticsSampler->sampleTemperature(&system);
        }
        // Allowing the system to equilibrate.
        system.setThermostatOn(false);
        for (int i = 0; i < 500; i++) {
            movie->saveState(&system);
            statisticsSampler->sample(&system, i);
            system.temperature = statisticsSampler->temperature();
            system.step(dt);
            statisticsSampler->sampleKineticEnergy(&system);
            statisticsSampler->sampleTemperature(&system);
        }
        system.save(filename);
    } else {
        // Loading the state and running the program.
        system.load(filename);
        for (int i = 0; i < timesteps; i++) {
            movie->saveState(&system);
            statisticsSampler->sample(&system, i);
            system.temperature = statisticsSampler->temperature();
            statisticsSampler->sampleKineticEnergy(&system);
            statisticsSampler->sampleTemperature(&system);
        }
        statisticsSampler->sampleHeatCapacity(&system);
    }

    movie->saveState(&system);
    statisticsSampler->closeFiles();
    movie->close();

    return 0;
}
