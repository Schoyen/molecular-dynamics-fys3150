from MDFramework import MDFramework

program = MDFramework(1e-15, 10, 1000, 100, 0.01, 1000, 500, 1000, 5.26, 2.5 * 3.405, 0)
program.run_MD_simulation()

# Allow MDFramework to decide if the program is compiled.
#program.compile_MD()
# Run this after the creation of article.
#program.clean_MD()
