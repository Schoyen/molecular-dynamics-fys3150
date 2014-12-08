from PlotPressureNumber import PlotPressureDensity
from os import walk

walk_list = [i for i in walk("doc/pressure")]
list_of_filenames1 = []
list_of_filenames2 = []

for i in range(len(walk_list[0][2])):
    if walk_list[0][2][i].startswith("pressure"):
        list_of_filenames1.append(walk_list[0][0] + "/" + walk_list[0][2][i])
    elif walk_list[0][2][i].startswith("number"):
        list_of_filenames2.append(walk_list[0][0] + "/" + walk_list[0][2][i])

list_of_pressures100 = []
list_of_pressures200 = []
list_of_pressures300 = []
for i in range(len(list_of_filenames1)):
    if list_of_filenames1[i].endswith("-100.txt"):
        list_of_pressures100.append(list_of_filenames1[i])
    elif list_of_filenames1[i].endswith("-200.txt"):
        list_of_pressures200.append(list_of_filenames1[i])
    else:
        list_of_pressures300.append(list_of_filenames1[i])

list_of_densities100 = []
list_of_densities200 = []
list_of_densities300 = []
for i in range(len(list_of_filenames2)):
    if list_of_filenames2[i].endswith("-100.txt"):
        list_of_densities100.append(list_of_filenames2[i])
    elif list_of_filenames2[i].endswith("-200.txt"):
        list_of_densities200.append(list_of_filenames2[i])
    else:
        list_of_densities300.append(list_of_filenames2[i])

p1 = PlotPressureDensity(list_of_pressures100, list_of_densities100)
p2 = PlotPressureDensity(list_of_pressures200, list_of_densities200)
p3 = PlotPressureDensity(list_of_pressures300, list_of_densities300)

p1.read_values()
p2.read_values()
p3.read_values()
p1.plot_values("Pressure - density for T = 100 K", "doc/PD100.pdf")
p2.plot_values("Pressure - density for T = 200 K", "doc/PD200.pdf")
p3.plot_values("Pressure - density for T = 300 K", "doc/PD300.pdf")
