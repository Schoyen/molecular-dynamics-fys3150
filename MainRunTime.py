from PlotRunTime import PlotRunTime
from os import walk

walk_list = [i for i in walk("doc/runtime")]
list_of_filenames1 = []
list_of_filenames2 = []

for i in range(len(walk_list[0][2])):
    if walk_list[0][2][i].endswith("-0.txt"):
        list_of_filenames1.append(walk_list[0][0] + "/" + walk_list[0][2][i])
    else:
        list_of_filenames2.append(walk_list[0][0] + "/" + walk_list[0][2][i])

p1 = PlotRunTime(list_of_filenames1)
p2 = PlotRunTime(list_of_filenames2)
p1.read_values()
p2.read_values()
p1.plot_values("Kilo atoms per timestep for cell lists", "doc/cellLists.pdf")
p2.plot_values("Kilo atoms per timestep for every atom pair", "doc/oldForce.pdf")
