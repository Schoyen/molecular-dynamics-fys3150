from matplotlib.pylab import plot, show, xlabel, ylabel, title, savefig

class PlotPressureDensity:

    def __init__(self, list_of_pressures_files, list_of_densities_files):
        self.list_of_pressures_files = list_of_pressures_files
        self.list_of_densities_files = list_of_densities_files

    def read_values(self):
        sum_pressures = 0
        self.list_of_pressures = []
        self.list_of_densities = []

        for i in self.list_of_pressures_files:
            with open(i, 'r') as f:
                counter = 0
                for line in f:
                    sum_pressures += float(line)
                    counter += 1
                self.list_of_pressures.append(sum_pressures / float(counter))
                counter = 0

        for i in self.list_of_densities_files:
            with open(i, 'r') as f:
                for line in f:
                    self.list_of_densities.append(float(line))

        self.list_of_densities = sorted(self.list_of_densities)
        self.list_of_densities.reverse()
        self.list_of_pressures = sorted(self.list_of_pressures)
        self.list_of_pressures.reverse()

    def plot_values(self, TITLE, SAVE):
        plot(self.list_of_densities, self.list_of_pressures)
        title(TITLE)
        xlabel("Densities")
        ylabel("Pressure")
        savefig(SAVE)
        show()
