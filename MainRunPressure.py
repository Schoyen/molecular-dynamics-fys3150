from MDFramework import MDFramework

program = MDFramework(1e-15, 10, 0, 0, 10, 1000, 1000, 0, 2.5 * 3.405, 0)
program.compile_MD()
program.measure_pressure()
