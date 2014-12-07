#include <iostream>
#include <fstream>
#include "math/random.h"
#include "potentials/potential.h"
#include "potentials/lennardjones.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include "berendsen.h"
#include <chrono>

using namespace std;
using namespace chrono;

int main()
{
    double dt = UnitConverter::timeFromSI(1e-14); // You should try different values for dt as well.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;

    auto start = high_resolution_clock::now();
    System system;
    // For more than 2 x 2 x 2 FCCLattice we need a bigger system size.
    int numberOfFCCLattices = 10;
    double cellSize = 2.5 * 3.405; // rcut
    int numberOfAtoms = 4 * numberOfFCCLattices * numberOfFCCLattices * numberOfFCCLattices;
    double initialTemperature = 100.0; // In Kelvin.
    system.createFCCLattice(numberOfFCCLattices, UnitConverter::lengthFromAngstroms(5.26), UnitConverter::temperatureFromSI(initialTemperature), cellSize);
    double tbath = 1000;
    double relaxationTime = 0.01; // Figure this one out.
    StatisticsSampler *statisticsSampler = new StatisticsSampler();
    system.setPotential(new LennardJones(3.405, 1.0));
    system.setIntegrator(new VelocityVerlet());
    system.setThermostat(new BerendsenThermostat(UnitConverter::temperatureFromSI(tbath), relaxationTime));
    system.removeMomentum();
    system.setForceCalculation(false);
    system.setThermostatOn(false);

    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    statisticsSampler->createFiles();
    string filename;
    //filename = "test.txt";
    //system.load(filename);
    for(int timestep=0; timestep<1000; timestep++) {
        movie->saveState(&system);
        system.temperature = statisticsSampler->temperature();
        system.step(dt);
        /*
        if (timestep < 300) {
            system.step(dt);
        } else {
            system.setThermostatOn(true); // Figure out something smarter.
            system.step(dt);
        }
        */
        statisticsSampler->sample(&system, timestep);
        statisticsSampler->sampleKineticEnergySquared(&system);
        statisticsSampler->sampleTotalKineticEnergy(&system);

        //std::cout << UnitConverter::energyToEv(statisticsSampler->totalEnergy()) << std::endl;

        cout << timestep << endl;
    }
    statisticsSampler->sampleHeatCapacity(&system);
    movie->saveState(&system);

    auto finish = high_resolution_clock::now();
    filename = "build/DATA/calculationTime" + to_string(numberOfAtoms) + ".txt";
    ofstream file(filename);
    if (file.is_open()) {
        file << duration_cast<nanoseconds>(finish - start).count()*(1.0e-9) << " s";
    } else cout << "Unable to write to file." << endl;


    file.close();
    statisticsSampler->closeFiles();
    movie->close();

    return 0;
}
