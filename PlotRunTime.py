from matplotlib.pylab import plot, show, title, xlabel, ylabel, savefig

class PlotRunTime:

    def __init__(self, list_of_filenames):
        # list_of_filenames should only contain one force calculation scheme.
        self.list_of_filenames = list_of_filenames

    def read_values(self):
        self.list_of_KAT = []
        self.list_of_time = []
        self.number_of_atoms = []

        index_dictionary = {}
        for counter, i in enumerate(self.list_of_filenames):
            with open(i, 'r') as f:
                for line in f:
                    self.list_of_KAT.append(float(line.split()[0]))
                    self.list_of_time.append(int(line.split()[1]))
                    self.number_of_atoms.append(int(line.split()[2]))
                    index_dictionary[self.number_of_atoms[counter]] = counter

        self.number_of_atoms = sorted(self.number_of_atoms)
        temp_list_of_KAT = []
        temp_list_of_time = []
        for i in self.number_of_atoms:
            temp_list_of_KAT.append(self.list_of_KAT[index_dictionary[i]])
            temp_list_of_time.append(self.list_of_time[index_dictionary[i]])

        self.list_of_KAT = temp_list_of_KAT
        self.list_of_time = temp_list_of_time


    def plot_values(self, TITLE, SAVE):
        plot(self.number_of_atoms, self.list_of_KAT)
        title(TITLE)
        xlabel("Number of atoms")
        ylabel("Kilo atoms per timestep")
        savefig(SAVE)
        show()

    def plot_time(self, TITLE, SAVE):
        plot(self.number_of_atoms, self.list_of_time)
        title(TITLE)
        xlabel("Number of atoms")
        ylabel("Time spent in seconds")
        savefig(SAVE)
        show()
