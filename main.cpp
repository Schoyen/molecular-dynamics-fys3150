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
    system.setSystemSize(UnitConverter::lengthFromAngstroms(vec3(30, 30, 30)));
    int numberOfFCCLattices = 5;
    double cellSize = 6;
    int numberOfAtoms = 4 * numberOfFCCLattices * numberOfFCCLattices * numberOfFCCLattices;
    system.createFCCLattice(numberOfFCCLattices, UnitConverter::lengthFromAngstroms(5.26), UnitConverter::temperatureFromSI(3000.0), cellSize);
    system.setPotential(new LennardJones(3.405, 1.0)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.removeMomentum();

    /*
    for(int n=0; n<100; n++) {
        // Add one example atom. You'll have to create many such atoms in the createFCCLattice function above.
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); // Argon mass, see http://en.wikipedia.org/wiki/Argon
        atom->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
        atom->position = vec3(0.0, 0.0, 0.0);
        //atom->position.randomUniform(0, system.systemSize().x());
        system.atoms().push_back(atom); // Add it to the list of atoms
    }
    */
    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //
    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    string filename;
    for(int timestep=0; timestep<1000; timestep++) {
        filename = "build/DATA/statistics" + to_string(timestep) + ".txt";
        system.step(dt);
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
