CXX = g++
RM = rm -r
CPPFLAGS = -Wno-write-strings -g -O3 -std=c++11 -pedantic -Wall
DEPS = atom.h io.h statisticssampler.h system.h unitconverter.h potentials/lennardjones.h potentials/potential.h math/random.h math/vec3.h integrators/eulercromer.h integrators/integrator.h integrators/velocityverlet.h celllist.h cell.h berendsen.h
TESTFLAGS = -lgtest
COMP = build/atom.o build/io.o build/statisticssampler.o build/system.o build/unitconverter.o build/POTENTIALS/lennardjones.o build/POTENTIALS/potential.o build/MATH/random.o build/MATH/vec3.o build/INTEGRATORS/eulercromer.o build/INTEGRATORS/integrator.o build/INTEGRATORS/velocityverlet.o build/celllist.o build/cell.o build/berendsen.o

build:
	mkdir -p build/MATH build/INTEGRATORS build/POTENTIALS build/DATA

build/atom.o: atom.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/io.o: io.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/statisticssampler.o: statisticssampler.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/system.o: system.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/unitconverter.o: unitconverter.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/celllist.o: celllist.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/cell.o: cell.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/berendsen.o: berendsen.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/POTENTIALS/lennardjones.o: potentials/lennardjones.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/POTENTIALS/potential.o: potentials/potential.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/MATH/random.o: math/random.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/MATH/vec3.o: math/vec3.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/INTEGRATORS/eulercromer.o: integrators/eulercromer.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/INTEGRATORS/integrator.o: integrators/integrator.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/INTEGRATORS/velocityverlet.o: integrators/velocityverlet.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $< -c -o $@

build/OLDMAINCPP: oldmain.cpp $(DEPS) $(COMP) | build
	$(CXX) $(CPPFLAGS) $(COMP) $< -o $@

clean:
	$(RM) build movie.xyz
