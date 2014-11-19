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
    double dt = UnitConverter::timeFromSI(1e-15); // You should try different values for dt as well.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;

    auto start = high_resolution_clock::now();
    System system;
    // For more than 2 x 2 x 2 FCCLattice we need a bigger system size.
    system.setSystemSize(UnitConverter::lengthFromAngstroms(vec3(28, 28, 28)));
    int numberOfFCCLattices = 5;
    double cellSize = 7;
    int numberOfAtoms = 4 * numberOfFCCLattices * numberOfFCCLattices * numberOfFCCLattices;
    double initialTemperature = 100.0; // In Kelvin.
    system.createFCCLattice(numberOfFCCLattices, UnitConverter::lengthFromAngstroms(5.26), UnitConverter::temperatureFromSI(initialTemperature), cellSize);
    double tbath = 3000;
    double relaxationTime = 0.01; // Figure this one out.
    system.setPotential(new LennardJones(3.405, 1.0, new BerendsenThermostat(UnitConverter::temperatureFromSI(tbath), relaxationTime))); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.removeMomentum();

    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //
    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    string filename;
    //filename = "test.txt";
    //system.load(filename);
    for(int timestep=0; timestep<200; timestep++) {
        if (timestep < 100) {
            filename = "build/DATA/statisticsTHERMO" + to_string(timestep) + ".txt";
            system.step(dt, true);
        } else {
            filename = "build/DATA/statistics" + to_string(timestep) + ".txt";
            system.step(dt, true);
            //system.save("test.txt");
            //break;
        }
        statisticsSampler->sample(&system, filename);

        movie->saveState(&system);
        //cout << timestep << endl;
    }
    auto finish = high_resolution_clock::now();
    filename = "build/DATA/calculationTime" + to_string(numberOfAtoms) + ".txt";
    ofstream file(filename);
    if (file.is_open()) {
        file << duration_cast<nanoseconds>(finish - start).count()*(1.0e-9) << " s";
    } else cout << "Unable to write to file." << endl;


    file.close();
    movie->close();

    return 0;
}
