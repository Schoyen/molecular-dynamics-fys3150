from os import path, system

class MDFramework:

    def __init__(self, dt, fcc_lattice, initial_temp, t_bath,\
                 relaxation_time, timesteps, timesteps_before_thermo, timesteps_after_thermo,\
                 lattice_constant, r_cut, old_force_calculation):
        self.dt = dt
        self.fcc_lattice = fcc_lattice
        self.initial_temp = initial_temp
        self.t_bath = t_bath
        self.relaxation_time = relaxation_time
        self.timesteps = timesteps

        # If the thermostat is going to be off all the time, set timesteps_before_thermo to -1.
        self.timesteps_before_thermo = timesteps_before_thermo
        self.timesteps_after_thermo = timesteps_after_thermo
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

    def create_article(self):
        # Use a script to put the documents in the right locations.
        # Let all files needed be created and let pdflatex run its course before cleaning.
        print ("""
==============================
Creating article and removing build.
==============================
""")

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
timesteps_before_thermo = %g
timesteps_after_thermo  = %g
lattice_constant        = %g
r_cut                   = %g
old_force_calculation   = %g
==============================
               """ % (self.dt, self.fcc_lattice, self.initial_temp,\
                      self.t_bath, self.relaxation_time, self.timesteps, self.timesteps_before_thermo,\
                      self.timesteps_after_thermo, self.lattice_constant, self.r_cut,\
                      self.old_force_calculation)) # *cringe*
        system("./build/MAINCPP %g %g %g %g %g %g %g %g %g %g %g" % (self.dt,\
               self.fcc_lattice, self.initial_temp, self.t_bath, self.relaxation_time,\
               self.timesteps, self.timesteps_before_thermo, self.timesteps_after_thermo,\
               self.lattice_constant, self.r_cut, self.old_force_calculation))

    def time_MD_simulation(self):
        # This method should time the two different force calculation methods and store the values as files.
        # Run the program for 1000 steps with 500 atoms.
        # Run the build/MAINCPP for different values.
        print("""
==============================
Timing program for cell lists
and calculation of forces
between all atom pairs.

We will run the program for
500 and 1000 atoms with 
100 - 1000 timesteps.

Time will be measured in 
kiloatom per timestep.

The temperature is irrelevant
during these calculations.
We will not use the
thermostat either.
==============================
              """)
        for i in range(0, 2, 1):
            for j in range(5, 11, 5):
                for k in range(100, 1001, 100):
                    system("./build/MAINCPP %g %g %g %g %g %g %g %g %g %g %g" % (self.dt,\
                           j, self.initial_temp, self.t_bath, self.relaxation_time,\
                           k, 0, 0,\
                           self.lattice_constant, self.r_cut, i))
                    system("sh sortRunTime.sh")
                    print("Done calculating %s for %g atoms with %g timesteps" % ("cell lists" if i == 0 else "force between all atom pairs",\
                                                                                  4 * j**3, k))
