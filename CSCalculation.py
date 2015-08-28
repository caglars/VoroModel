__author__ = 'caglars'

import random
import tess
import math
import numpy
import pandas as pd
from scipy import stats

import CSSolution
import CSPlot
import CSProperties
import CSData

class CSCalculator():
    def __init__(self):
        self.betaConstant = 1.127
        self.alphaConstant = 5.615

        #self.particles = 25
        #self.x_min = 0
        #self.x_max = 6000
        #self.y_min = 0
        #self.y_max = 6001
        #self.z_min = 0
        #self.z_max = 60
        #self.boxLimits = [[self.x_min, self.y_min, self.z_min],
        #                  [self.x_max, self.y_max, self.z_max]]
        #self.pressure = []
        #self.pressure = [6000 for x in range(self.particles)]
        #self.permeabilityX = [0.01 for x in range(self.particles)]
        #self.permeabilityY = [0.01 for x in range(self.particles)]
        #self.permeabilityZ = [0.01 for x in range(self.particles)]
        #self.porosity = [0.18 for x in range(self.particles)]
        #self.referencePressure = 3031.
        #self.referenceDepth = 7000. # the top most depth
        pass

    def initialize(self, aContainer):
        referencePressure = self.reader.readSingleFloatValue("REFPRES")
        self.pressure = [referencePressure for x in range(self.particles)]
        self.finalPressure = [referencePressure for x in range(self.particles)]
        myProperties = CSProperties.CSFluidProperties()
        # TODO the calculation of the initial pressure should be iterative
        refDepth = myProperties.referenceDepth
        refPressure = myProperties.referencePressure

        newGammaFluid = 0.
        for cell in aContainer:
            #print("cell %s" % cell.id)
            gammaFluid = myProperties.findGammaFluid(self.pressure[cell.id])
            while True:
                depthOfCell = refDepth + aContainer[cell.id].pos[2]
                self.finalPressure[cell.id] = refPressure + gammaFluid*(depthOfCell - refDepth)
                newGammaFluid = myProperties.findGammaFluid(self.finalPressure[cell.id])
                #print("gamma %s and newGamma %s" % (gammaFluid, newGammaFluid))
                if abs(newGammaFluid - gammaFluid) > 0.001:
                    gammaFluid = newGammaFluid
                else:
                    break

            pass

        for x in range(0, self.particles):
            self.pressure[x] = self.finalPressure[x]-0.1
            pass


        #print(self.pressure)

    def rnd(self, myBoxLimits):
        ''' 1 is added and subtracted to prevent an error caused by random chooses exactly the limit values '''

        point_list = [(random.randint(myBoxLimits[0][0] + 1, myBoxLimits[1][0] - 1)),
                     (random.randint(myBoxLimits[0][1] + 1, myBoxLimits[1][1] - 1)),
                     (random.randint(myBoxLimits[0][2] + 1, myBoxLimits[1][2] - 1))]
        return point_list

    def readData(self):

        self.permeabilityX = self.reader.readValues("PERMX")
        self.permeabilityY = self.reader.readValues("PERMY")
        self.permeabilityZ = self.reader.readValues("PERMZ")
        self.porosity = self.reader.readValues("PORO")
        self.wellRadius = self.reader.readWellData(dataType="RW")
        self.wellSkins = self.reader.readWellData(dataType="SKIN")
        self.productionRate = self.reader.readWellData(dataType = "FLOWRATE")
        self.wellThickness = self.reader.readWellData(dataType = "PERFTHICK")



    def buildModel(self):

        self.reader = CSData.CSDataReader()


        cellList =self.reader.readParticles()
        self.particles = len(cellList)

        containerLimits = self.reader.readValues("LIMITS")
        self.x_min = 0
        self.x_max = containerLimits[0]
        self.y_min = 0
        self.y_max = containerLimits[1]
        self.z_min = 0
        self.z_max = containerLimits[2]


        cntr = tess.Container(cellList, limits=[(self.x_min, self.y_min, self.z_min),
                                                (self.x_max, self.y_max, self.z_max)], periodic=False)

        return cntr


    def findDistance(self, coordParticle1, coordParticle2):
        return math.sqrt(math.pow((coordParticle2[0] - coordParticle1[0]), 2)
                     + math.pow((coordParticle2[1] - coordParticle1[1]), 2)
                     + math.pow((coordParticle2[2] - coordParticle1[2]), 2))



    def getPermeability(self, aDataFrame, cellId, neighborId):
        cellPosition = aDataFrame['positions'][cellId]
        neighborPosition = aDataFrame['positions'][neighborId]
        cellVolume = aDataFrame['volumes'][cellId]
        neighborVolume = aDataFrame['volumes'][neighborId]

        deltaX = abs(cellPosition[0] - neighborPosition[0])
        deltaY = abs(cellPosition[1] - neighborPosition[1])
        deltaZ = abs(cellPosition[2] - neighborPosition[2])


        if deltaX == 0:
            if deltaY == 0:
                if deltaZ == 0:
                    #cell and neighbor same (not likely to occur this)
                    #print("cell and neighbor are same!!!")
                    return 0
                # z only
                #print ("deltaX and deltaY = 0")
                permAverageZ = (cellVolume*self.permeabilityZ[cellId]+neighborVolume*self.permeabilityZ[neighborId])/(cellVolume+neighborVolume)
                return permAverageZ
            if deltaZ == 0:
                # y only
                #print ("deltaX and deltaZ = 0")
                permAverageY = (cellVolume*self.permeabilityY[cellId]+neighborVolume*self.permeabilityY[neighborId])/(cellVolume+neighborVolume)
                return permAverageY
            # y and z
            radian = math.atan(deltaZ/deltaY)
            permAverageY = (cellVolume*self.permeabilityY[cellId]+neighborVolume*self.permeabilityY[neighborId])/(cellVolume+neighborVolume)
            permAverageZ = (cellVolume*self.permeabilityZ[cellId]+neighborVolume*self.permeabilityZ[neighborId])/(cellVolume+neighborVolume)
            #print("deltaX = 0")
            #print("cell %s and neigh %s" % (cellId, neighborId))
            #print("cell volume %s and neighbor volume %s" % (cellVolume, neighborVolume))
            #print("cell permY %s and neighbor permY %s" % (self.permeabilityY[cellId],self.permeabilityY[neighborId]))
            #print("cell permZ %s and neighbor permZ %s" % (self.permeabilityZ[cellId],self.permeabilityZ[neighborId]))
            #print("avgPermY %s and avgPermZ %s" % (permAverageY, permAverageZ))
            #print ("radian %s" % radian)
            permYPrime = math.cos(radian)*permAverageY
            permZPrime = math.sin(radian)*permAverageZ
            #permYPrime = stats.hmean([math.cos(radian)*self.permeabilityY[cellId], math.cos(radian)*self.permeabilityY[neighborId]])
            #permZPrime = stats.hmean([math.sin(radian)*self.permeabilityZ[cellId], math.sin(radian)*self.permeabilityZ[neighborId]])
            #print("yP %s zP %s" % (permYPrime, permZPrime))
            return math.sqrt(permYPrime**2 + permZPrime**2)

        elif deltaY == 0:
            if deltaZ == 0:
                # x only
                #print ("deltaY and deltaZ = 0")
                permAverageX = (cellVolume*self.permeabilityX[cellId]+neighborVolume*self.permeabilityX[neighborId])/(cellVolume+neighborVolume)
                return permAverageX
            # x and z
            radian = math.atan(deltaZ/deltaX)
            permAverageX = (cellVolume*self.permeabilityX[cellId]+neighborVolume*self.permeabilityX[neighborId])/(cellVolume+neighborVolume)
            permAverageZ = (cellVolume*self.permeabilityZ[cellId]+neighborVolume*self.permeabilityZ[neighborId])/(cellVolume+neighborVolume)
            #print("deltaY = 0")
            #print("cell %s and neigh %s" % (cellId, neighborId))
            #print("cell volume %s and neighbor volume %s" % (cellVolume, neighborVolume))
            #print("cell permX %s and neighbor permX %s" % (self.permeabilityX[cellId],self.permeabilityX[neighborId]))
            #print("cell permZ %s and neighbor permZ %s" % (self.permeabilityZ[cellId],self.permeabilityZ[neighborId]))
            #print("avgPermX %s and avgPermZ %s" % (permAverageX, permAverageZ))
            #print ("radian %s" % radian)
            permXPrime = math.cos(radian)*permAverageX
            permZPrime = math.sin(radian)*permAverageZ
            #permYPrime = stats.hmean([math.cos(radian)*self.permeabilityY[cellId], math.cos(radian)*self.permeabilityY[neighborId]])
            #permZPrime = stats.hmean([math.sin(radian)*self.permeabilityZ[cellId], math.sin(radian)*self.permeabilityZ[neighborId]])
            #print("xP %s zP %s" % (permXPrime, permZPrime))
            return math.sqrt(permXPrime**2 + permZPrime**2)

        elif deltaZ == 0:
            # x and y
            radian = math.atan(deltaY/deltaX)
            permAverageX = (cellVolume*self.permeabilityX[cellId]+neighborVolume*self.permeabilityX[neighborId])/(cellVolume+neighborVolume)
            permAverageY = (cellVolume*self.permeabilityY[cellId]+neighborVolume*self.permeabilityY[neighborId])/(cellVolume+neighborVolume)
            #print("deltaZ = 0")
            #print("cell %s and neigh %s" % (cellId, neighborId))
            #print("cell volume %s and neighbor volume %s" % (cellVolume, neighborVolume))
            #print("cell permX %s and neighbor permX %s" % (self.permeabilityX[cellId],self.permeabilityX[neighborId]))
            #print("cell permY %s and neighbor permY %s" % (self.permeabilityY[cellId],self.permeabilityY[neighborId]))
            #print("avgPermX %s and avgPermY %s" % (permAverageX, permAverageY))
            #print ("radian %s" % radian)
            permXPrime = math.cos(radian)*permAverageX
            permYPrime = math.sin(radian)*permAverageY
            #permYPrime = stats.hmean([math.cos(radian)*self.permeabilityY[cellId], math.cos(radian)*self.permeabilityY[neighborId]])
            #permZPrime = stats.hmean([math.sin(radian)*self.permeabilityZ[cellId], math.sin(radian)*self.permeabilityZ[neighborId]])
            #print("xP %s yP %s" % (permXPrime, permYPrime))
            return math.sqrt(permXPrime**2 + permYPrime**2)
        else:
            #deltaX, deltaY and deltaZ are not zero or all of them are zero. All zero should not come to this function
            #xy and z
            #first we need to find xy
            radianXY = math.atan(deltaY/deltaX)
            permAverageX = (cellVolume*self.permeabilityX[cellId]+neighborVolume*self.permeabilityX[neighborId])/(cellVolume+neighborVolume)
            permAverageY = (cellVolume*self.permeabilityY[cellId]+neighborVolume*self.permeabilityY[neighborId])/(cellVolume+neighborVolume)
            permAverageZ = (cellVolume*self.permeabilityZ[cellId]+neighborVolume*self.permeabilityZ[neighborId])/(cellVolume+neighborVolume)
            #print("Nothing is zero")
            #print("cell %s and neigh %s" % (cellId, neighborId))
            #print("cell volume %s and neighbor volume %s" % (cellVolume, neighborVolume))
            #print("cell permX %s and neighbor permX %s" % (self.permeabilityX[cellId],self.permeabilityX[neighborId]))
            #print("cell permY %s and neighbor permY %s" % (self.permeabilityY[cellId],self.permeabilityY[neighborId]))
            #print("avgPermX %s and avgPermY %s" % (permAverageX, permAverageY))
            #print("deltaX %s deltaY %s delta Z %s" % (deltaX, deltaY, deltaZ))
            #print ("radianXY %s" % radianXY)
            permXPrime = math.cos(radianXY)*permAverageX
            permYPrime = math.sin(radianXY)*permAverageY
            #permYPrime = stats.hmean([math.cos(radian)*self.permeabilityY[cellId], math.cos(radian)*self.permeabilityY[neighborId]])
            #permZPrime = stats.hmean([math.sin(radian)*self.permeabilityZ[cellId], math.sin(radian)*self.permeabilityZ[neighborId]])
            #print("xP %s yP %s" % (permXPrime, permYPrime))
            permXY = math.sqrt(permXPrime**2 + permYPrime**2)

            lengthXY = math.sqrt(deltaX**2 + deltaY**2)
            radian = math.atan(deltaZ/lengthXY)
            permZPrime = math.sin(radian)*permAverageZ
            permXYPrime = math.cos(radian)*permXY
            #print("cell permZ %s and neighbor permZ %s" % (self.permeabilityZ[cellId],self.permeabilityZ[neighborId]))
            #print("permXY %s" % permXY)
            #print("radian %s" % radian)
            #print("xyP %s zP %s" % (permXYPrime, permZPrime))
            return math.sqrt(permZPrime**2 + permXYPrime**2)


        #print("x %s y %s z %s" % (deltaX, deltaY, deltaZ))




    def getFormationVolumeFactor(self, aDataFrame, cellId, neighborId):
        # explicit treatment of the parameters, pressures are at time step n
        cellVolume = aDataFrame['volumes'][cellId]
        neighborVolume = aDataFrame['volumes'][neighborId]
        #todo Factor is calculated using bulk volumes, it would be better to use pore volumes
        factor = cellVolume/(cellVolume+neighborVolume)
        myProperty = CSProperties.CSFluidProperties()
        formationVolumeFactor = factor*myProperty.findFormationVolumeFactor(self.finalPressure[cellId]) \
                                + (1-factor)*myProperty.findFormationVolumeFactor(self.finalPressure[neighborId])
        return formationVolumeFactor
        #return 1.0

    def getViscosity(self, aDataFrame, cellId, neighborId):
        # explicit treatment of the parameters, pressures are at time step n
        cellVolume = aDataFrame['volumes'][cellId]
        neighborVolume = aDataFrame['volumes'][neighborId]
        #todo Factor is calculated using bulk volumes, it would be better to use pore volumes
        factor = cellVolume/(cellVolume+neighborVolume)
        myProperty = CSProperties.CSFluidProperties()
        viscosity = factor*myProperty.findViscosity(self.finalPressure[cellId],1) \
                                + (1-factor)*myProperty.findViscosity(self.finalPressure[neighborId],1)
        return viscosity
        #return 10.

    # TODO It would be better to read all the data at the beginning and check for errors instead of reading in separate functions
    def getBottomHolePressure(self, aContainer, cellId, blockPressure):
        # Using Eqn 6.32 of Ertekin
        #self.wellRadius = numpy.zeros(self.particles)
        #self.wellSkins = numpy.zeros(self.particles)
        myProperty = CSProperties.CSFluidProperties()
        viscosity = myProperty.findViscosity(blockPressure, 1)
        formationVolumeFactor = myProperty.findFormationVolumeFactor(blockPressure)
        equivalentRadius = math.pow((3*aContainer[cellId].volume()/(4*math.pi)), 1./3.)
        harmonicPermeability = math.pow((self.permeabilityX[cellId]*self.permeabilityY[cellId]*self.permeabilityZ[cellId]),1./3.)
        #print("viscosity: %s" % viscosity)
        #print("blockPressure: %s" % blockPressure)
        #print("FVF: %s" % formationVolumeFactor)
        #print("req: %s" % equivalentRadius)
        #print("kH: %s" % harmonicPermeability)
        #print("rw: %s" % self.wellRadius[cellId])
        #print("s: %s" % self.wellSkins[cellId])
        #print("pi: %s" % math.pi)
        #print("h: %s" % self.wellThickness[cellId])

        bottomHolePressure = blockPressure + self.productionRate[cellId]*viscosity*formationVolumeFactor\
                                             *(math.log(equivalentRadius/self.wellRadius[cellId]+self.wellSkins[cellId]-0.5))\
                                             /(2*math.pi*self.betaConstant*harmonicPermeability*self.wellThickness[cellId])
        return bottomHolePressure





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


    def simRunSlightlyCompressible(self, aContainer):
        #self.particles = numberOfParticles

        mySolver = CSSolution.CSSolver()
        myProperties = CSProperties.CSFluidProperties()
        myPlotter = CSPlot.CSPlotter(self.particles, self.pressure)

        #viscosity = self.reader.readSingleFloatValue("VISCOSITY")
        #formationVolumeFactor = self.reader.readSingleFloatValue("FORMATIONVOLUMEFACTOR")
        #liquidCompressibility = self.reader.readSingleFloatValue("LIQUIDCOMPRESSIBILITY")
        #referenceFormationVolumeFactor = self.reader.readSingleFloatValue("REFFVF")
        deltaTime = self.reader.readSingleFloatValue("DELTATIME")
        numberOfTimeSteps = self.reader.readSingleIntValue("TIMESTEPS")

        #viscosity = 10
        #formationVolumeFactor = 1
        #liquidCompressibility = 3.5E-6
        #referenceFormationVolumeFactor = 1
        #fluidDensity = 62.4
        #gravityAcceleration = 32.17

        #alphaConstant = 5.615
        #betaConstant = 1.127
        #gammaConstant = 0.21584e-3

        #deltaTime = 15
        #length = 0
        #totalCoefficient = 0

        #numberOfTimeSteps = 1

        x = 0
        y = 0
        z = 0

        #print(self.reader.readWellData(dataType = "FLOWRATE"))

        #self.productionRate = numpy.zeros(self.particles)


        gamma = numpy.zeros(self.particles)
        gravity = numpy.zeros((self.particles, self.particles))
        rightHandSide = numpy.zeros(self.particles)
        coefficient = numpy.zeros((self.particles, self.particles))

        referencePressure = self.reader.readSingleFloatValue("REFPRES")
        gammaFluid = myProperties.findGammaFluid(referencePressure)
        #print (gammaFluid)

        #productionRate[particles - 1] = -150.0
        #productionRate[particles - 2] = -200.0
        #productionRate[self.particles - 1] = -100.0

        #self.productionRate = self.reader.readWellData(dataType = "FLOWRATE")

        neighborSeries = pd.Series()
        posSeries = pd.Series()
        faceAreaSeries = pd.Series()
        volumeSeries = pd.Series()

        for cell in aContainer:
            neighborSeries = neighborSeries.append(pd.Series([cell.neighbors()], index=[cell.id]))
            posSeries = posSeries.append(pd.Series([cell.pos], index=[cell.id]))
            faceAreaSeries = faceAreaSeries.append(pd.Series([cell.face_areas()], index=[cell.id]))
            volumeSeries = volumeSeries.append(pd.Series([cell.volume()], index=[cell.id]))

        data = {'neighbors': neighborSeries,
                'positions': posSeries,
                'faceAreas': faceAreaSeries,
                'volumes': volumeSeries}
        myDataFrame = pd.DataFrame(data)

        for timeStep in range(0, numberOfTimeSteps):
            print("timeStep: %s" % timeStep)
            for iteration in range(0, 1):
                print("iteration: %s" % iteration)
                for cellIndex in range(0, len(myDataFrame)):
                    if cellIndex%100 == 0:
                        print("cellIndex: %s" % cellIndex)
                    totalCoefficient = 0
                    totalGravity = 0
                    totalRightHandSideGravity = 0
                    numberOfNeighbors = len(myDataFrame['neighbors'][cellIndex])
                    for neighborIndex in range(0, numberOfNeighbors):
                        neighbor = myDataFrame['neighbors'][cellIndex][neighborIndex]
                        if neighbor >= 0:
                            length = self.findDistance(myDataFrame['positions'][cellIndex],myDataFrame['positions'][neighbor])
                            self.permeability = self.getPermeability(myDataFrame, cellIndex, neighbor)
                            coefficient[cellIndex][neighbor] = (self.betaConstant
                                                                * myDataFrame['faceAreas'][cellIndex][neighborIndex]
                                                                * self.permeability) / (self.getViscosity(myDataFrame, cellIndex, neighbor)
                                                                  * self.getFormationVolumeFactor(myDataFrame, cellIndex, neighbor)
                                                                  * length)


                            # TODO gammaFluid should change for each neighbor according to their pressure and density
                            gravity[cellIndex][neighbor] = gammaFluid * coefficient[cellIndex][neighbor]
                            totalRightHandSideGravity += gravity[cellIndex][neighbor]*myDataFrame['positions'][neighbor][2]
                            totalGravity += gravity[cellIndex][neighbor]
                            totalCoefficient += coefficient[cellIndex][neighbor]
                            pass
                        pass
                    # Since the fluid is slightly compressible,
                    # TODO the gamma value is calculated using equation 8.133 in Ertekin but porosity change should be included
                    # this version assumes constant porosity

                    if (self.pressure[cellIndex]-self.finalPressure[cellIndex]) == 0.0:
                        #print("cellID: %s - yes denominator is zero" % cell.id)
                        self.pressure[cellIndex]=self.finalPressure[cellIndex]-0.1

                    gamma[cellIndex] = (myDataFrame['volumes'][cellIndex]/self.alphaConstant) \
                                     * ((self.porosity[cellIndex]/myProperties.findFormationVolumeFactor(self.finalPressure[cellIndex]))
                                        * (myProperties.findFormationVolumeFactor(self.finalPressure[cellIndex])
                                           /myProperties.findFormationVolumeFactor(self.pressure[cellIndex])-1)
                                        /(self.pressure[cellIndex]-self.finalPressure[cellIndex]))

                    coefficient[cellIndex][cellIndex] = -totalCoefficient - (gamma[cellIndex] / deltaTime)

                    # TODO The gravity CG value should be checked from the equation
                    gravity[cellIndex][cellIndex] = -totalGravity

                    totalRightHandSideGravity += gravity[cellIndex][cellIndex]*myDataFrame['positions'][cellIndex][2]
                    #print("total RHS gravity %s" % totalRightHandSideGravity)
                    rightHandSide[cellIndex] = -(self.productionRate[cellIndex]
                                               + (gamma[cellIndex] / deltaTime) * self.finalPressure[cellIndex]
                                               - totalRightHandSideGravity)

                    pass

                print("sending to numpy")
                solutionArray = mySolver.numpySolver(coefficient, rightHandSide)

                for x in range(0, self.particles):
                    self.pressure[x] = solutionArray[x]
                    pass

                myPlotter.pressureAtParticleDetail(3, timeStep, iteration, self.pressure)

            print("Time step: %s" % timeStep)
            print("pressure is written to pressureOut.csv")

            for x in range(0, self.particles):
                self.finalPressure[x] = self.pressure[x]
                pass

            #silinecek sonra
            bottomHolePressure = self.getBottomHolePressure(aContainer, 3, self.finalPressure[x])
            myPlotter.writeBottomHolePressure(x,timeStep,bottomHolePressure)


            myPlotter.gnuplot(aContainer)

            myPlotter.graphr(aContainer, timeStep)

            myPlotter.pressureAtParticle(3, timeStep, self.pressure)

            #myPlotter.pressureAtParticle(24, timeStep, self.pressure)

            myPlotter.permeabilityGraphr(aContainer, self.permeabilityX, self.permeabilityY, self.permeabilityZ)


            pass




        # todo An error approach is required to end the iterations
        for timeStep in range(0, numberOfTimeSteps):
            #print("timestep %s" % timeStep)
            '''
            if timeStep == 29:
                deltaTime = 2
                self.productionRate[3] = 0.
            '''
            for iteration in range(0, 1):
                #print("iteration %s" % iteration)
                for cell in aContainer:
                    # print("cellid: %s, pressure: %s finalPressure: %s" % (cell.id, self.pressure[cell.id], self.finalPressure[cell.id]))
                    #print(cell.id)
                    neighborCounter = 0
                    totalCoefficient = 0
                    totalGravity = 0
                    totalRightHandSideGravity = 0
                    for neighbor in cell.neighbors():
                        # print("neighbor %s and face_area = %s" % (neighbor, cell.face_areas()[neighborCounter]))
                        if neighbor >= 0:
                            length = self.findDistance(aContainer[cell.id].pos, aContainer[neighbor].pos)
                            self.permeability = self.getPermeability(aContainer, cell.id, neighbor)
                            #print (self.permeability)
                            coefficient[cell.id][neighbor] = (self.betaConstant * cell.face_areas()[neighborCounter]
                                                              * self.permeability
                                                              / (self.getViscosity(aContainer, cell.id, neighbor)
                                                                 * self.getFormationVolumeFactor(aContainer, cell.id, neighbor) * length))
                            #print(coefficient[cell.id][neighbor])
                            #print("area %s" % cell.face_areas()[neighborCounter])

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
                    # TODO the gamma value is calculated using equation 8.133 in Ertekin but porosity change should be included
                    # this version assumes constant porosity


                    if (self.pressure[cell.id]-self.finalPressure[cell.id]) == 0.0:
                        #print("cellid: %s - yes denominator is zero" % cell.id)
                        self.pressure[cell.id]=self.finalPressure[cell.id]-0.1

                    gamma[cell.id] = (cell.volume()/self.alphaConstant) \
                                     * ((self.porosity[cell.id]/myProperties.findFormationVolumeFactor(self.finalPressure[cell.id]))
                                        * (myProperties.findFormationVolumeFactor(self.finalPressure[cell.id])
                                           /myProperties.findFormationVolumeFactor(self.pressure[cell.id])-1)
                                        /(self.pressure[cell.id]-self.finalPressure[cell.id]))

                    '''
                    if cell.id==578:
                        #print("vol: %s FVFold: %s FVFnew: %s pres: %s finPres: %s" % (cell.volume(), myProperties.findFormationVolumeFactor(self.finalPressure[cell.id]), myProperties.findFormationVolumeFactor(self.pressure[cell.id]),self.pressure[cell.id], self.finalPressure[cell.id]))
                        print("time step: %s iteration: %s pressure: %s" % (timeStep, iteration, self.pressure[cell.id]))
                    '''


                    coefficient[cell.id][cell.id] = -totalCoefficient - (gamma[cell.id] / deltaTime)
                    #print(coefficient[cell.id][cell.id])

                    # TODO The gravity CG value should be checked from the equation
                    gravity[cell.id][cell.id] = -totalGravity

                    totalRightHandSideGravity += gravity[cell.id][cell.id]*aContainer[cell.id].pos[2]
                    #print("total RHS gravity %s" % totalRightHandSideGravity)
                    rightHandSide[cell.id] = -(self.productionRate[cell.id]
                                               + (gamma[cell.id] / deltaTime) * self.finalPressure[cell.id]
                                               - totalRightHandSideGravity)
                    #print(rightHandSide)

                    # print("cell id = %s volume = %s" % (cell.id, cell.volume()))
                    pass

                # solutionArray = mySolver.simpleSolver(particles, coefficient, rightHandSide, pressure)

                '''
                tempFile = open("coefficient{}.dat".format(iteration), "w")
                temp2File = open("righthandside{}.dat".format(iteration), "w")


                print(len(coefficient))
                print(len(rightHandSide))

                for a in range(0, self.particles):
                    for b in range(0, self.particles):
                        if coefficient[a][b] != 0.0:
                            tempFile.write("a: %s b:%s %s \n" % (a, b, str(coefficient[a][b])))

                for a in range (0, self.particles):
                    temp2File.write("a: %s %s \n" % (a, str(rightHandSide[a])))

                tempFile.close()
                temp2File.close()
                '''




                solutionArray = mySolver.numpySolver(coefficient, rightHandSide)

                for x in range(0, self.particles):
                    self.pressure[x] = solutionArray[x]
                    pass

                myPlotter.pressureAtParticleDetail(3, timeStep, iteration, self.pressure)

            print("Time step: %s" % timeStep)
            print("pressure is written to pressureOut.csv")

            '''
            for x in range(0, self.particles):
                self.finalPressure[x] = self.pressure[x] - 0.1
                pass
            '''

            for x in range(0, self.particles):
                self.finalPressure[x] = self.pressure[x]
                pass

            #active kuyular için shut-in olsa dahi gerekebilir
            '''
            for x in range(0, self.particles):
                if self.productionRate[x] != 0:
                    bottomHolePressure = self.getBottomHolePressure(aContainer, x, self.finalPressure[x])
                    myPlotter.writeBottomHolePressure(x,timeStep,bottomHolePressure)
            '''

            #silinecek sonra
            bottomHolePressure = self.getBottomHolePressure(aContainer, 3, self.finalPressure[x])
            myPlotter.writeBottomHolePressure(x,timeStep,bottomHolePressure)


            myPlotter.gnuplot(aContainer)

            myPlotter.graphr(aContainer, timeStep)

            myPlotter.pressureAtParticle(3, timeStep, self.pressure)

            #myPlotter.pressureAtParticle(24, timeStep, self.pressure)

            myPlotter.permeabilityGraphr(aContainer, self.permeabilityX, self.permeabilityY, self.permeabilityZ)


            pass