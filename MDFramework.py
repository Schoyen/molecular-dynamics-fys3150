from os import path, system
from matplotlib.pylab import plot, show, xlabel, ylabel, savefig, title

class MDFramework:

    def __init__(self, dt, fcc_lattice, initial_temp, t_bath,\
                 relaxation_time, timesteps, timesteps_with_thermo,\
                 lattice_constant, r_cut, old_force_calculation):
        self.dt = dt
        self.fcc_lattice = fcc_lattice
        self.initial_temp = initial_temp
        self.t_bath = t_bath
        self.relaxation_time = relaxation_time
        self.timesteps = timesteps

        self.timesteps_with_thermo = timesteps_with_thermo
        self.lattice_constant = lattice_constant
        self.r_cut = r_cut
        self.old_force_calculation = old_force_calculation

    def compile_MD(self):
        if not path.exists("build"):
            print ("""
==============================
Compiling program.
==============================
""")
            system("make")
            system("make build/MAINCPP")

    def clean_MD(self):
        print ("""
==============================
Removing build.
==============================
""")
        if path.exists("build"):
            system("make clean")

    def run_MD_simulation(self):
        print ("""
==============================
Running program for:
dt                      = %g
fcc_lattice             = %g
initial_temp            = %g
t_bath                  = %g
relaxation_time         = %g
timesteps               = %g
timesteps_with_thermo   = %g
lattice_constant        = %g
r_cut                   = %g
old_force_calculation   = %g
==============================
               """ % (self.dt, self.fcc_lattice, self.initial_temp,\
                      self.t_bath, self.relaxation_time, self.timesteps,\
                      self.timesteps_with_thermo, self.lattice_constant, self.r_cut,\
                      self.old_force_calculation)) # *cringe*
        system("./build/MAINCPP %g %g %g %g %g %g %g %g %g %g" % (self.dt,\
               self.fcc_lattice, self.initial_temp, self.t_bath, self.relaxation_time,\
               self.timesteps, self.timesteps_with_thermo,\
               self.lattice_constant, self.r_cut, self.old_force_calculation))

    def initial_behaviour(self):
        self.run_MD_simulation()
        """
        system("mv build/DATA/time.txt doc/initialbehaviour/")
        system("mv build/DATA/kineticEnergy.txt doc/initialbehaviour/")
        system("mv build/DATA/potentialEnergy.txt doc/initialbehaviour/")
        system("mv build/DATA/totalEnergy.txt doc/initialbehaviour/")
        system("mv build/DATA/pressure.txt doc/initialbehaviour/")
        system("mv build/DATA/temperature.txt doc/initialbehaviour/")
        """
        
        time = self.reader("doc/initialbehaviour/time.txt")
        kinetic = self.reader("doc/initialbehaviour/kineticEnergy.txt")
        potential = self.reader("doc/initialbehaviour/potentialEnergy.txt")
        total = self.reader("doc/initialbehaviour/totalEnergy.txt")
        pressure = self.reader("doc/initialbehaviour/pressure.txt")
        temperature = self.reader("doc/initialbehaviour/temperature.txt")
        n = len(time)

        self.plotter(time, kinetic, "time [s]", "kinetic energy [eV]", "Kinetic energy as a function of time", "doc/kineticEnergy.pdf")
        self.plotter(time, potential, "time [s]", "potential energy [eV]", "Potential energy as a function of time", "doc/potentialEnergy.pdf")
        self.plotter(time, total, "time [s]", "total energy [eV]", "Total energy as a function of time", "doc/totalEnergy.pdf")
        self.plotter(time, pressure, "time [s]", "pressure [N/m^2]", "Pressure as a function of time", "doc/pressure.pdf")
        self.plotter(time, temperature, "time [s]", "temperature [K]", "Temperature as a function of time", "doc/temperature.pdf")

    def reader(self, filename):
        tempList = []
        with open(filename, 'r') as f:
            for line in f:
                tempList.append(float(line))

        return tempList

    def plotter(self, first_val, second_val, XLABEL, YLABEL,  TITLE, SAVE):
        plot(first_val, second_val)
        xlabel(XLABEL)
        ylabel(YLABEL)
        title(TITLE)
        savefig(SAVE)
        show()


    def measure_pressure(self):
        # Measure for different temperatures and lattice constants.
        # Save the heat capacity as well. Check against experimental values.
        print("""
==============================
Measuring the pressure of 
the system for different
values of the lattice 
constant b and different 
temperatures.

Using the thermostat to
achieve the desired 
temperature and then 
allowing the system to 
equilibrate.
==============================
              """)
        b_list = [4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
        for b in b_list:
            for T in range(100, 301, 100):
                self.initial_temp = T
                self.t_bath = T
                self.lattice_constant = b
                self.use_thermostat()
                system("mv build/DATA/pressure.txt doc/pressure/pressure-" + str(b) + "-" + str(T) + ".txt")
                system("mv build/DATA/numberDensity.txt doc/pressure/numberDensity-" + str(b) + "-" + str(T) + ".txt")
                system("mv build/DATA/heatCapacity.txt doc/pressure/heatCapacity-" + str(b) + "-" + str(T) + ".txt")
            print("Done calculating for b = %g" % b)

    def use_thermostat(self):
        # Run the program with a thermostat then save it when it is done.
        # Load the new saved program and then run without thermostat.

        # The last 0 or 1 tells us if the thermostat should be off or on respectively.
        system("./build/MAINCPP %g %g %g %g %g %g %g %g %g %g %g" % (self.dt,\
               self.fcc_lattice, self.initial_temp, self.t_bath, self.relaxation_time,\
               self.timesteps, self.timesteps_with_thermo,\
               self.lattice_constant, self.r_cut, self.old_force_calculation, 1))
        system("./build/MAINCPP %g %g %g %g %g %g %g %g %g %g %g" % (self.dt,\
               self.fcc_lattice, self.initial_temp, self.t_bath, self.relaxation_time,\
               self.timesteps, self.timesteps_with_thermo,\
               self.lattice_constant, self.r_cut, self.old_force_calculation, 0))

    def time_MD_simulation(self, limit):
        print("""
==============================
Timing program for cell lists
and calculation of forces
between all atom pairs.

We will run the program for
1000 timesteps and
varying number of atoms.

Time will be measured in 
kiloatom per timestep.

The temperature is irrelevant
during these calculations.
We will not use the
thermostat either.
==============================
              """)
        for i in range(0, 2, 1):
            for j in range(6, limit, 1):
                system("./build/MAINCPP %g %g %g %g %g %g %g %g %g %g" % (self.dt,\
                       j, self.initial_temp, self.t_bath, self.relaxation_time,\
                       1000, 0,\
                       self.lattice_constant, self.r_cut, i))
                system("sh sortRunTime.sh")
                print("Done calculating %s for %g atoms with %g timesteps" % ("force with cell lists" if i == 0 else "force between all atom pairs",\
                                                                              4 * j**3, 1000))
