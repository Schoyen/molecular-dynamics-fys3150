#include <iostream>
#include <cstdlib>
#include "io.h"
#include "system.h"
#include "unitconverter.h"

using namespace std;
/*
 * List of things that need to be sent to the main-method.
 * double dt
 * int systemSize (only quadratic cubes)
 * int numberOfFCCLattices
 * int cellSize
 * double initialTemperature
 * double tbath
 * double relaxationTime
 * int timestep
 * int timestepStartThermostat
 * int timestepEndThermostat
 * double latticeConstant
 * bool oldForceCalculation
 */
int main(int argc, char* argv[])
{
    if (argc < 13) {
        cout << "\n==================================" << endl;
        cout << "Not enough command line arguments." << endl;
        cout << "Usage: " << argv[0] << " 1 2 3 4 5 6 7 8 9 10 11 12\n"
             << "1: double dt\n"
             << "2: int system size (only quadratic cubes)\n"
             << "3: int number of FCC lattices\n"
             << "4: int cell size\n"
             << "5: double initial temperature\n"
             << "6: double t_bath\n"
             << "7: double relaxation time\n"
             << "8: int number of timesteps\n"
             << "9: int number of timesteps before turning on the thermostat\n"
             << "10: int number of timesteps before turning off the thermostat\n"
             << "11: bool 1='true' for old force calculation and 0 (or anything) ='false' for cell lists"
             << endl;
        cout << "==================================\n" << endl;
        exit(1);
    }

    double dt = atof(argv[1]);
    int systemSize = atoi(argv[2]);
    int numberOfFCCLattices = atoi(argv[3]);
    int cellSize = atoi(argv[4]);
    double initialTemperature = atof(argv[5]);
    double tbath = atof(argv[6]);
    double relaxationTime = atof(argv[7]);
    int timestep = atoi(argv[8]);
    int timestepStartThermostat = atoi(argv[9]);
    int timestepEndThermostat = atoi(argv[10]);
    bool oldForceCalculation;
    if (atoi(argv[11]) == 1) oldForceCalculation = true;
    else oldForceCalculation = false;

    IO *movie = new IO();
    System system;
    system.setSystemSize(UnitConverter::lengthFromAngstroms(vec3(systemSize, systemSize, systemSize)));
    system.createFCCLattice(numberOfFCCLattices, UnitConverter::lengthFromAngstroms(

    if (oldeForceCalculation) {
    } else {
    }
}
