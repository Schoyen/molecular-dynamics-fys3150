from MDFramework import MDFramework

program = MDFramework(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
program.run_MD_simulation()

# Allow MDFramework to decide if the program is compiled.
#program.compile_MD()
# Run this after the creation of article.
#program.clean_MD()
