from numpy import zeros
from matplotlib.pylab import title, plot, show, xlabel, ylabel, hold, legend, savefig

class ForceCalculation:

    def __init__(self, sigma, epsilon, distance):
        self.sigma = sigma
        self.epsilon = epsilon
        self.distance = distance
        self.n = self.distance.size

    def calculateForce(self):
        self.force = zeros(self.n)
        for i in range(self.n):
            self.force[i] = self.distance[i] * (4 * self.epsilon * ((self.sigma**12/float(self.distance[i]**14)) - (self.sigma**6/float(self.distance[i]**8))))

    def plotForce(self, TITLE):
        plot(self.distance, self.force)
        hold('on')
        shorterDistance = zeros(25)
        for i in range(25):
            shorterDistance[i] = self.distance[4*i]
        plot(shorterDistance, zeros(25), '.')
        legend(('force', 'zero axis',), loc=1)
        xlabel("relative distance r_ij")
        ylabel("force F(r_ij)")
        title(TITLE)
        savefig('doc/forcePlot.png')
        show()
