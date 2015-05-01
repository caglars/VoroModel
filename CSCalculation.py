__author__ = 'caglars'

import random
import tess
import math

import CSSolution
import CSPlot

class CSCalculator():
    def __init__(self):
        self.particles = 50
        self.x_min = -1
        self.x_max = 6001
        self.y_min = -1
        self.y_max = 6001
        self.z_min = -1
        self.z_max = 81
        self.boxLimits = [[self.x_min, self.y_min, self.z_min],
                          [self.x_max, self.y_max, self.z_max]]
        self.pressure = [6000 for x in range(self.particles)]
        pass

    def rnd(self, myBoxLimits):
        ''' 1 is added and subtracted to prevent an error caused by random chooses exactly the limit values '''

        point_list = [(random.randint(myBoxLimits[0][0] + 1, myBoxLimits[1][0] - 1)),
                     (random.randint(myBoxLimits[0][1] + 1, myBoxLimits[1][1] - 1)),
                     (random.randint(myBoxLimits[0][2] + 1, myBoxLimits[1][2] - 1))]
        return point_list

    def readData(self):

        self.permeability = self.readDataFor('permeability.dat')

        print("After getting permeability")

        print(self.permeability)

        self.pressure = self.readDataFor('initialPressure.dat')

        print("After getting pressure")

        print (self.pressure)



    def readDataFor(self, fileName):
        self.temp = []
        self.valueList = []

        with open(fileName, 'r') as f:
            data = f.read()

        f.close()

        words = data.split()

        print(words)

        for text in words:
            self.temp.append(text.split('*'))

        print(self.temp)

        for count, value in enumerate(self.temp):
            print(count)
            print(value)
            if len(value) == 1:
                print(float(value[0]))
                self.valueList.append(float(value[0]))
            else:
                print('no %s' % float(value[1]))
                for count in range(0, int(value[0])):
                    self.valueList.append(float(value[1]))

        print(self.valueList)

        self.temp.clear()

        return self.valueList

    def calculate(self, numberOfParticles):

        self.particles = numberOfParticles
        # i = None
        # x, y, z, r = None
        cellList = []

        '''
        cellList.append([1000, 1000, 40])
        cellList.append([1000, 2000, 40])
        cellList.append([1000, 3000, 40])
        cellList.append([1000, 4000, 40])
        cellList.append([1000, 5000, 40])
        cellList.append([2000, 1000, 40])
        cellList.append([2000, 2000, 40])
        cellList.append([2000, 3000, 40])
        cellList.append([2000, 4000, 40])
        cellList.append([2000, 5000, 40])
        cellList.append([3000, 1000, 40])
        cellList.append([3000, 2000, 40])
        cellList.append([3000, 4000, 40])
        cellList.append([3000, 5000, 40])
        cellList.append([4000, 1000, 40])
        cellList.append([4000, 2000, 40])
        cellList.append([4000, 3000, 40])
        cellList.append([4000, 4000, 40])
        cellList.append([4000, 5000, 40])
        cellList.append([5000, 1000, 40])
        cellList.append([5000, 2000, 40])
        cellList.append([5000, 3000, 40])
        cellList.append([5000, 4000, 40])
        cellList.append([5000, 5000, 40])
        cellList.append([3000, 3000, 40])

        '''

        for x in range(0, self.particles-1):
            cellList.append(self.rnd(self.boxLimits))
            pass
        cellList.append([3000, 3000, 40])



        cntr = tess.Container(cellList, limits=[(self.x_min, self.y_min, self.z_min),
                                                (self.x_max, self.y_max, self.z_max)],
                              periodic=False)

        return cntr


    def distance(self, coordParticle1, coordParticle2):
        return math.sqrt(math.pow((coordParticle2[0] - coordParticle1[0]), 2)
                     + math.pow((coordParticle2[1] - coordParticle1[1]), 2)
                     + math.pow((coordParticle2[2] - coordParticle1[2]), 2))


    def simRun(self, aContainer, numberOfParticles):
        self.particles = numberOfParticles

        mySolver = CSSolution.CSSolver()

        # permeability = 0.015
        viscosity = 10
        formationVolumeFactor = 1
        porosity = 0.18
        liquidCompressibility = 3.5E-6
        referansFormationVolumeFactor = 1
        alphaConstant = 5.615
        betaConstant = 1.127
        deltaTime = 15
        length = 0
        totalCoefficient = 0

        numberOfTimeSteps = 2

        x = 0
        y = 0
        z = 0

        productionRate = [0 for x in range(self.particles)]
        gamma = [0 for x in range(self.particles)]
        # self.pressure = [6000 for x in range(self.particles)]
        rightHandSide = [0 for x in range(self.particles)]
        coefficient = [[0 for x in range(self.particles)] for x in range(self.particles)]

        #productionRate[particles - 1] = -150.0
        #productionRate[particles - 2] = -200.0
        productionRate[self.particles - 1] = -100.0

        for timeStep in range(0, numberOfTimeSteps):
            for cell in aContainer:
                neighborCounter = 0
                totalCoefficient = 0
                for neighbor in cell.neighbors():
                    # print("neighbor %s and face_area = %s" % (neighbor, cell.face_areas()[neighborCounter]))
                    if neighbor >= 0:
                        length = self.distance(aContainer[cell.id].pos, aContainer[neighbor].pos)
                        coefficient[cell.id][neighbor] = betaConstant * cell.face_areas()[
                            neighborCounter] * self.permeability[cell.id] / (viscosity * formationVolumeFactor * length)
                        totalCoefficient = totalCoefficient + coefficient[cell.id][neighbor]
                        # print("i = %s, neighbor = %s, neighborCounter = %s" % (cell.id, neighbor, neighborCounter))

                        pass
                    pass
                    neighborCounter = neighborCounter + 1

                gamma[cell.id] = (cell.volume() * porosity * liquidCompressibility) / (
                    alphaConstant * referansFormationVolumeFactor)
                coefficient[cell.id][cell.id] = -totalCoefficient - gamma[cell.id] / deltaTime

                rightHandSide[cell.id] = -productionRate[cell.id] - (gamma[cell.id] / deltaTime) * self.pressure[cell.id]

                # print("cell id = %s volume = %s" % (cell.id, cell.volume()))
                pass

            # solutionArray = mySolver.simpleSolver(particles, coefficient, rightHandSide, pressure)
            solutionArray = mySolver.numpySolver(coefficient, rightHandSide)

            for x in range(0, self.particles):
                self.pressure[x] = solutionArray[x]
                pass

            print("Time step: %s" % timeStep)
            print("pressure is written to pressureOut.csv")

            myPlotter = CSPlot.CSPlotter(self.particles, self.pressure)

            myPlotter.gnuplot(aContainer)

            myPlotter.graphr(aContainer)

            return self.pressure

            pass

