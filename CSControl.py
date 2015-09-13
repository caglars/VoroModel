__author__ = 'caglars'


#Porosity change is ignored
#Gravity is added
#pressure dependant FVF and viscosity are added
#Working for undersaturated reservoir slightly compressible fluid


import time

import CSCalculation

class CSController:
    def __init__(self):
        self.particles = 25



    def main(self):
        #self.pressure = [6000 for x in range(self.particles)]
        self.respond()

    def respond(self):

        t1 = time.time()

        calculator = CSCalculation.CSCalculator()

        cont = calculator.buildModel()

        calculator.readData()

        calculator.initialize(cont)

        #calculator.simRunIncompressible(cont, self.particles)

        calculator.simRunSlightlyCompressible(cont)

        t2 = time.time()
        print('it took %s seconds' % (t2 - t1))

run = CSController()
run.main()