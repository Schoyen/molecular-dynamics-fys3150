from MDFramework import MDFramework

program = MDFramework(1e-15, 10, 300, 1000, 1, 1000, 1000, 5.26, 2.5 * 3.405, 1)
program.compile_MD()
program.time_MD_simulation(31)
