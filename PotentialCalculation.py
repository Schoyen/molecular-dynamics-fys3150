from numpy import zeros
from matplotlib.pylab import title, plot, show, xlabel, ylabel, hold, legend, savefig

class PotentialCalculation:

    def __init__(self, sigma, epsilon, distance):
        self.sigma = sigma
        self.epsilon = epsilon
        self.distance = distance
        self.n = self.distance.size

    def calculatePotential(self):
        self.potential = zeros(self.n)
        for i in range(self.n):
            self.potential[i] = 4 * self.epsilon * ((self.sigma / float(self.distance[i]))**(12) -\
                                                    (self.sigma / float(self.distance[i]))**6)

    def plotPotential(self, TITLE):
        plot(self.distance/self.sigma, self.potential/self.epsilon)
        hold('on')
        shorterDistance = zeros(25)
        for i in range(25):
            shorterDistance[i] = self.distance[4*i]
        plot(shorterDistance/self.sigma, zeros(25), '.')
        legend(('potential', 'zero axis',), loc=1)
        xlabel("relative distance r_ij/sigma")
        ylabel("potential U(r_ij)/epsilon")
        title(TITLE)
        savefig('doc/potentialPlot.png')
        show()
