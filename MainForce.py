from ForceCalculation import ForceCalculation
from numpy import zeros

distance = zeros(100)
sigma = 3.405
epsilon = 1.0

for i in range(100):
    distance[i] = 0.1*i + (5.26/2.0)

fc = ForceCalculation(sigma, epsilon, distance)
fc.calculateForce()
fc.plotForce("Reduction of force as a function of relative distance.")
