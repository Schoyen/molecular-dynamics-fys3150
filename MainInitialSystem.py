from MDFramework import MDFramework

program = MDFramework(1e-15, 10, 300, 1, 1, 2000, 0, 5.26, 2.5 * 3.405, 0)
program.compile_MD()
program.initial_behaviour()
