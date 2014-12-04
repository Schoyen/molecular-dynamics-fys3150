from PotentialCalculation import PotentialCalculation
from numpy import zeros

distance = zeros(100)
sigma = 3.405
epsilon = 1.0

for i in range(5, 105):
    distance[i-5] = 0.1*i + 5.26/2.0

pc = PotentialCalculation(sigma, epsilon, distance)
pc.calculatePotential()
pc.plotPotential("Reduction of potential as a function of relative distance.")
