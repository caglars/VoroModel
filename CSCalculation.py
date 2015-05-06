__author__ = 'caglars'

import random
import tess
import math
import numpy
from scipy import stats

import CSSolution
import CSPlot

class CSCalculator():
    def __init__(self):
        self.particles = 25
        self.x_min = -1
        self.x_max = 6001
        self.y_min = -1
        self.y_max = 6001
        self.z_min = -1
        self.z_max = 61
        self.boxLimits = [[self.x_min, self.y_min, self.z_min],
                          [self.x_max, self.y_max, self.z_max]]
        self.pressure = [6000 for x in range(self.particles)]
        self.permeabilityX = [0.01 for x in range(self.particles)]
        self.permeabilityY = [0.01 for x in range(self.particles)]
        self.permeabilityZ = [0.01 for x in range(self.particles)]
        self.porosity = [0.18 for x in range(self.particles)]
        pass

    def rnd(self, myBoxLimits):
        ''' 1 is added and subtracted to prevent an error caused by random chooses exactly the limit values '''

        point_list = [(random.randint(myBoxLimits[0][0] + 1, myBoxLimits[1][0] - 1)),
                     (random.randint(myBoxLimits[0][1] + 1, myBoxLimits[1][1] - 1)),
                     (random.randint(myBoxLimits[0][2] + 1, myBoxLimits[1][2] - 1))]
        return point_list

    def readData(self):

        self.permeabilityX = self.readDataFor('permeabilityX.dat')
        self.permeabilityY = self.readDataFor('permeabilityY.dat')
        self.permeabilityZ = self.readDataFor('permeabilityZ.dat')
        self.pressure = self.readDataFor('initialPressure.dat')
        self.porosity = self.readDataFor('porosity.dat')



    def readDataFor(self, fileName):
        self.temp = []
        self.valueList = []

        with open(fileName, 'r') as f:
            data = f.read()

        f.close()

        words = data.split()

        #print(words)

        for text in words:
            self.temp.append(text.split('*'))

        #print(self.temp)

        for count, value in enumerate(self.temp):
            #print(count)
            #print(value)
            if len(value) == 1:
                #print(float(value[0]))
                self.valueList.append(float(value[0]))
            else:
                #print('no %s' % float(value[1]))
                for count in range(0, int(value[0])):
                    self.valueList.append(float(value[1]))

        #print(self.valueList)

        self.temp.clear()

        return self.valueList

    def calculate(self, numberOfParticles):

        self.particles = numberOfParticles
        # i = None
        # x, y, z, r = None
        cellList = []


        cellList.append([1000, 1000, 50])
        cellList.append([1000, 2000, 50])
        cellList.append([1000, 3000, 50])
        cellList.append([1000, 4000, 50])
        cellList.append([1000, 5000, 50])
        cellList.append([2000, 1000, 40])
        cellList.append([2000, 2000, 40])
        cellList.append([2000, 3000, 40])
        cellList.append([2000, 4000, 40])
        cellList.append([2000, 5000, 40])
        cellList.append([3000, 1000, 30])
        cellList.append([3000, 2000, 30])
        cellList.append([3000, 4000, 30])
        cellList.append([3000, 5000, 30])
        cellList.append([4000, 1000, 20])
        cellList.append([4000, 2000, 20])
        cellList.append([4000, 3000, 20])
        cellList.append([4000, 4000, 20])
        cellList.append([4000, 5000, 20])
        cellList.append([5000, 1000, 10])
        cellList.append([5000, 2000, 10])
        cellList.append([5000, 3000, 10])
        cellList.append([5000, 4000, 10])
        cellList.append([5000, 5000, 10])
        cellList.append([3000, 3000, 30])

        '''

        for x in range(0, self.particles-1):
            cellList.append(self.rnd(self.boxLimits))
            pass
        cellList.append([3000, 3000, 40])


        cellList.append([1000, 5000, 10])
        cellList.append([2000, 4000, 20])
        cellList.append([3000, 3000, 20])
        cellList.append([4000, 2000, 40])
        cellList.append([5000, 1000, 50])

        '''


        cntr = tess.Container(cellList, limits=[(self.x_min, self.y_min, self.z_min),
                                                (self.x_max, self.y_max, self.z_max)],
                              periodic=False)

        return cntr


    def distance(self, coordParticle1, coordParticle2):
        return math.sqrt(math.pow((coordParticle2[0] - coordParticle1[0]), 2)
                     + math.pow((coordParticle2[1] - coordParticle1[1]), 2)
                     + math.pow((coordParticle2[2] - coordParticle1[2]), 2))

    def getPermeability(self, aContainer, cellId, neighborId):
        cellPosition = aContainer[cellId].pos
        neighborPosition = aContainer[neighborId].pos
        cellVolume = aContainer[cellId].volume()
        neighborVolume = aContainer[neighborId].volume()
        deltaX = abs(cellPosition[0] - neighborPosition[0])
        deltaY = abs(cellPosition[1] - neighborPosition[1])
        deltaZ = abs(cellPosition[2] - neighborPosition[2])


        if deltaX == 0:
            if deltaY == 0:
                if deltaZ == 0:
                    #cell and neighbor same (not likely to occur this)
                    print("cell and neighbor are same!!!")
                    return 0
                # z only
                print ("deltaX and deltaY = 0")
                permAverageZ = (cellVolume*self.permeabilityZ[cellId]+neighborVolume*self.permeabilityZ[neighborId])/(cellVolume+neighborVolume)
                return permAverageZ
            if deltaZ == 0:
                # y only
                print ("deltaX and deltaZ = 0")
                permAverageY = (cellVolume*self.permeabilityY[cellId]+neighborVolume*self.permeabilityY[neighborId])/(cellVolume+neighborVolume)
                return permAverageY
            # y and z
            radian = math.atan(deltaZ/deltaY)
            permAverageY = (cellVolume*self.permeabilityY[cellId]+neighborVolume*self.permeabilityY[neighborId])/(cellVolume+neighborVolume)
            permAverageZ = (cellVolume*self.permeabilityZ[cellId]+neighborVolume*self.permeabilityZ[neighborId])/(cellVolume+neighborVolume)
            print("deltaX = 0")
            print("cell %s and neigh %s" % (cellId, neighborId))
            print("cell volume %s and neighbor volume %s" % (cellVolume, neighborVolume))
            print("cell permY %s and neighbor permY %s" % (self.permeabilityY[cellId],self.permeabilityY[neighborId]))
            print("cell permZ %s and neighbor permZ %s" % (self.permeabilityZ[cellId],self.permeabilityZ[neighborId]))
            print("avgPermY %s and avgPermZ %s" % (permAverageY, permAverageZ))
            print ("radian %s" % radian)
            permYPrime = math.cos(radian)*permAverageY
            permZPrime = math.sin(radian)*permAverageZ
            #permYPrime = stats.hmean([math.cos(radian)*self.permeabilityY[cellId], math.cos(radian)*self.permeabilityY[neighborId]])
            #permZPrime = stats.hmean([math.sin(radian)*self.permeabilityZ[cellId], math.sin(radian)*self.permeabilityZ[neighborId]])
            print("yP %s zP %s" % (permYPrime, permZPrime))
            return math.sqrt(permYPrime**2 + permZPrime**2)

        elif deltaY == 0:
            if deltaZ == 0:
                # x only
                print ("deltaY and deltaZ = 0")
                permAverageX = (cellVolume*self.permeabilityX[cellId]+neighborVolume*self.permeabilityX[neighborId])/(cellVolume+neighborVolume)
                return permAverageX
            # x and z
            radian = math.atan(deltaZ/deltaX)
            permAverageX = (cellVolume*self.permeabilityX[cellId]+neighborVolume*self.permeabilityX[neighborId])/(cellVolume+neighborVolume)
            permAverageZ = (cellVolume*self.permeabilityZ[cellId]+neighborVolume*self.permeabilityZ[neighborId])/(cellVolume+neighborVolume)
            print("deltaY = 0")
            print("cell %s and neigh %s" % (cellId, neighborId))
            print("cell volume %s and neighbor volume %s" % (cellVolume, neighborVolume))
            print("cell permX %s and neighbor permX %s" % (self.permeabilityX[cellId],self.permeabilityX[neighborId]))
            print("cell permZ %s and neighbor permZ %s" % (self.permeabilityZ[cellId],self.permeabilityZ[neighborId]))
            print("avgPermX %s and avgPermZ %s" % (permAverageX, permAverageZ))
            print ("radian %s" % radian)
            permXPrime = math.cos(radian)*permAverageX
            permZPrime = math.sin(radian)*permAverageZ
            #permYPrime = stats.hmean([math.cos(radian)*self.permeabilityY[cellId], math.cos(radian)*self.permeabilityY[neighborId]])
            #permZPrime = stats.hmean([math.sin(radian)*self.permeabilityZ[cellId], math.sin(radian)*self.permeabilityZ[neighborId]])
            print("xP %s zP %s" % (permXPrime, permZPrime))
            return math.sqrt(permXPrime**2 + permZPrime**2)

        elif deltaZ == 0:
            # x and y
            radian = math.atan(deltaY/deltaX)
            permAverageX = (cellVolume*self.permeabilityX[cellId]+neighborVolume*self.permeabilityX[neighborId])/(cellVolume+neighborVolume)
            permAverageY = (cellVolume*self.permeabilityY[cellId]+neighborVolume*self.permeabilityY[neighborId])/(cellVolume+neighborVolume)
            print("deltaZ = 0")
            print("cell %s and neigh %s" % (cellId, neighborId))
            print("cell volume %s and neighbor volume %s" % (cellVolume, neighborVolume))
            print("cell permX %s and neighbor permX %s" % (self.permeabilityX[cellId],self.permeabilityX[neighborId]))
            print("cell permY %s and neighbor permY %s" % (self.permeabilityY[cellId],self.permeabilityY[neighborId]))
            print("avgPermX %s and avgPermY %s" % (permAverageX, permAverageY))
            print ("radian %s" % radian)
            permXPrime = math.cos(radian)*permAverageX
            permYPrime = math.sin(radian)*permAverageY
            #permYPrime = stats.hmean([math.cos(radian)*self.permeabilityY[cellId], math.cos(radian)*self.permeabilityY[neighborId]])
            #permZPrime = stats.hmean([math.sin(radian)*self.permeabilityZ[cellId], math.sin(radian)*self.permeabilityZ[neighborId]])
            print("xP %s yP %s" % (permXPrime, permYPrime))
            return math.sqrt(permXPrime**2 + permYPrime**2)
        else:
            #deltaX, deltaY and deltaZ are not zero or all of them are zero. All zero should not come to this function
            #xy and z
            #first we need to find xy
            radianXY = math.atan(deltaY/deltaX)
            permAverageX = (cellVolume*self.permeabilityX[cellId]+neighborVolume*self.permeabilityX[neighborId])/(cellVolume+neighborVolume)
            permAverageY = (cellVolume*self.permeabilityY[cellId]+neighborVolume*self.permeabilityY[neighborId])/(cellVolume+neighborVolume)
            permAverageZ = (cellVolume*self.permeabilityZ[cellId]+neighborVolume*self.permeabilityZ[neighborId])/(cellVolume+neighborVolume)
            print("Nothing is zero")
            print("cell %s and neigh %s" % (cellId, neighborId))
            print("cell volume %s and neighbor volume %s" % (cellVolume, neighborVolume))
            print("cell permX %s and neighbor permX %s" % (self.permeabilityX[cellId],self.permeabilityX[neighborId]))
            print("cell permY %s and neighbor permY %s" % (self.permeabilityY[cellId],self.permeabilityY[neighborId]))
            print("avgPermX %s and avgPermY %s" % (permAverageX, permAverageY))
            print("deltaX %s deltaY %s delta Z %s" % (deltaX, deltaY, deltaZ))
            print ("radianXY %s" % radianXY)
            permXPrime = math.cos(radianXY)*permAverageX
            permYPrime = math.sin(radianXY)*permAverageY
            #permYPrime = stats.hmean([math.cos(radian)*self.permeabilityY[cellId], math.cos(radian)*self.permeabilityY[neighborId]])
            #permZPrime = stats.hmean([math.sin(radian)*self.permeabilityZ[cellId], math.sin(radian)*self.permeabilityZ[neighborId]])
            print("xP %s yP %s" % (permXPrime, permYPrime))
            permXY = math.sqrt(permXPrime**2 + permYPrime**2)

            lengthXY = math.sqrt(deltaX**2 + deltaY**2)
            radian = math.atan(deltaZ/lengthXY)
            permZPrime = math.sin(radian)*permAverageZ
            permXYPrime = math.cos(radian)*permXY
            print("cell permZ %s and neighbor permZ %s" % (self.permeabilityZ[cellId],self.permeabilityZ[neighborId]))
            print("permXY %s" % permXY)
            print("radian %s" % radian)
            print("xyP %s zP %s" % (permXYPrime, permZPrime))
            return math.sqrt(permZPrime**2 + permXYPrime**2)


        #print("x %s y %s z %s" % (deltaX, deltaY, deltaZ))
        return 0

    #TODO simRunIncompressible should be updated
    def simRunIncompressible(self, aContainer, numberOfParticles):
        self.particles = numberOfParticles

        mySolver = CSSolution.CSSolver()

        viscosity = 10
        formationVolumeFactor = 1
        liquidCompressibility = 3.5E-6
        referansFormationVolumeFactor = 1
        alphaConstant = 5.615
        betaConstant = 1.127
        deltaTime = 15
        length = 0
        totalCoefficient = 0

        numberOfTimeSteps = 10



        #productionRate = [0 for x in range(self.particles)]
        #gamma = [0 for x in range(self.particles)]
        #rightHandSide = [0 for x in range(self.particles)]
        #coefficient = [[0 for x in range(self.particles)] for x in range(self.particles)]
        productionRate = numpy.zeros(self.particles)
        gamma = numpy.zeros(self.particles)
        rightHandSide = numpy.zeros(self.particles)
        coefficient = numpy.zeros((self.particles, self.particles))

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
                        coefficient[cell.id][neighbor] = (betaConstant * cell.face_areas()[neighborCounter]
                                                          * self.permeabilityX[cell.id]
                                                          / (viscosity * formationVolumeFactor * length))
                        totalCoefficient = totalCoefficient + coefficient[cell.id][neighbor]
                        # print("i = %s, neighbor = %s, neighborCounter = %s" % (cell.id, neighbor, neighborCounter))

                        pass
                    pass
                    neighborCounter = neighborCounter + 1

                gamma[cell.id] = ((cell.volume() * self.porosity[cell.id] * liquidCompressibility)
                                  / (alphaConstant * referansFormationVolumeFactor))

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

            myPlotter.graphr(aContainer, timeStep)

            pass


    def simRunSlightlyCompressible(self, aContainer, numberOfParticles):
        self.particles = numberOfParticles

        mySolver = CSSolution.CSSolver()

        viscosity = 10
        formationVolumeFactor = 1
        liquidCompressibility = 3.5E-6
        referansFormationVolumeFactor = 1
        fluidDensity = 62.4
        gravityAcceleration = 32.17

        alphaConstant = 5.615
        betaConstant = 1.127
        gammaConstant = 0.21584e-3

        deltaTime = 15
        length = 0
        totalCoefficient = 0

        numberOfTimeSteps = 100

        x = 0
        y = 0
        z = 0


        productionRate = numpy.zeros(self.particles)
        gamma = numpy.zeros(self.particles)
        gravity = numpy.zeros((self.particles, self.particles))
        rightHandSide = numpy.zeros(self.particles)
        coefficient = numpy.zeros((self.particles, self.particles))


        # The density of the fluid should change and so this parameter
        # TODO I need to calculate this value as pressure and density changes
        gammaFluid = gammaConstant*fluidDensity*gravityAcceleration
        print (gammaFluid)

        #productionRate[particles - 1] = -150.0
        #productionRate[particles - 2] = -200.0
        productionRate[self.particles - 1] = -100.0

        for timeStep in range(0, numberOfTimeSteps):
            for cell in aContainer:
                neighborCounter = 0
                totalCoefficient = 0
                totalGravity = 0
                totalRightHandSideGravity = 0
                for neighbor in cell.neighbors():
                    # print("neighbor %s and face_area = %s" % (neighbor, cell.face_areas()[neighborCounter]))
                    if neighbor >= 0:
                        length = self.distance(aContainer[cell.id].pos, aContainer[neighbor].pos)
                        self.permeability = self.getPermeability(aContainer, cell.id, neighbor)
                        print (self.permeability)
                        coefficient[cell.id][neighbor] = (betaConstant * cell.face_areas()[neighborCounter]
                                                          * self.permeability
                                                          / (viscosity * formationVolumeFactor * length))

                        # TODO gammaFluid should change for each neighbor according to their pressure and density
                        gravity[cell.id][neighbor] = gammaFluid * coefficient[cell.id][neighbor]
                        totalRightHandSideGravity += gravity[cell.id][neighbor]*aContainer[neighbor].pos[2]
                        totalGravity += gravity[cell.id][neighbor]
                        totalCoefficient += coefficient[cell.id][neighbor]
                        # print("i = %s, neighbor = %s, neighborCounter = %s" % (cell.id, neighbor, neighborCounter))

                        pass
                    pass
                    neighborCounter = neighborCounter + 1


                # Since the fluid is slightly compressible,
                # TODO the gamma value should be calculated using equation 8.94 in Page 188
                # this version assumes constant porosity

                gamma[cell.id] = ((cell.volume() * self.porosity[cell.id] * liquidCompressibility)
                                  / (alphaConstant * referansFormationVolumeFactor))


                coefficient[cell.id][cell.id] = -totalCoefficient - (gamma[cell.id] / deltaTime)

                # TODO The gravity CG value should be checked from the equation
                gravity[cell.id][cell.id] = -totalGravity

                totalRightHandSideGravity += gravity[cell.id][cell.id]*aContainer[cell.id].pos[2]
                print("total RHS gravity %s" % totalRightHandSideGravity)
                rightHandSide[cell.id] = -(productionRate[cell.id]
                                           + (gamma[cell.id] / deltaTime) * self.pressure[cell.id]
                                           - totalRightHandSideGravity)

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

            myPlotter.graphr(aContainer, timeStep)

            myPlotter.pressureAtParticle(4, timeStep, self.pressure)

            myPlotter.permeabilityGraphr(aContainer, self.permeabilityX, self.permeabilityY, self.permeabilityZ)


            pass

