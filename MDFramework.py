from os import path, system

class MDFramework:

    def __init__(self, dt, system_size, fcc_lattice, cell_size, initial_temp, t_bath,\
                 relaxation_time, timesteps, timesteps_before_thermo, timesteps_after_thermo,\
                 old_force_calculation):
        self.dt = dt
        self.system_size = system_size
        self.fcc_lattice = fcc_lattice
        self.cell_size = cell_size
        self.initial_temp = initial_temp
        self.t_bath = t_bath
        self.relaxation_time = relaxation_time
        self.timesteps = timesteps
        self.timesteps_before_thermo = timesteps_before_thermo
        self.timesteps_after_thermo = timesteps_after_thermo
        self.old_force_calculation = old_force_calculation

    def compile_MD(self):
        if not path.exists("build"):
            print ("Compiling program.")
            system("make")
            system("make build/MAINCPP")

    def clean_MD(self):
        print ("Removing build.")
        if path.exists("build"):
            system("make clean")

    def create_article(self):
        # Let all files needed be created and let pdflatex run its course before cleaning.
        print ("Creating article and removing build.")

    def run_MD_simulation(self):
        if self.old_force_calculation == 1:
            self.time_MD_simulation()
        else:
            self.compile_MD()
            print ("""
==============================
Running program for:
dt                      = %g
system_size             = %g
fcc_lattice             = %g
cell_size               = %g
initial_temp            = %g
t_bath                  = %g
relaxation_time         = %g
timesteps               = %g
timesteps_before_thermo = %g
timesteps_after_thermo  = %g
old_force_calculation   = %g
==============================
                   """ % (self.dt, self.system_size, self.fcc_lattice, self.cell_size, self.initial_temp,\
                          self.t_bath, self.relaxation_time, self.timesteps, self.timesteps_before_thermo,\
                          self.timesteps_after_thermo, self.old_force_calculation)) # *cringe*
            system("./build/MAINCPP %g %g %g %g %g %g %g %g %g %g %g" % (self.dt, self.system_size,\
                   self.fcc_lattice, self.cell_size, self.initial_temp, self.t_bath, self.relaxation_time,\
                   self.timesteps, self.timesteps_before_thermo, self.timesteps_after_thermo,\
                   self.old_force_calculation))

    def time_MD_simulation(self):
        # This method should time the two different force calculation methods and store the values as files.
        # Run the program for 1000 steps with 500 atoms.
        pass
