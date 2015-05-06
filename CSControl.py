__author__ = 'caglars'

import time

import CSCalculation
import CSPlot

class CSController:
    def __init__(self):
        self.particles = 25



    def main(self):
        #self.pressure = [6000 for x in range(self.particles)]
        self.respond()

    def respond(self):

        t1 = time.time()

        calculator = CSCalculation.CSCalculator()

        cont = calculator.calculate(self.particles)

        calculator.readData()

        calculator.initialize(cont)

        #calculator.simRunIncompressible(cont, self.particles)

        calculator.simRunSlightlyCompressible(cont, self.particles)

        t2 = time.time()
        print('it took %s seconds' % (t2 - t1))

run = CSController()
run.main()