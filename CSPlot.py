__author__ = 'caglars'

import pandas as pd
import os

class CSPlotter():


    def __init__(self, numberOfParticles, pressure):
        #global ESCAPE
        #global numberOfParticles
        #global rotx, roty, beginx, beginy, rotate


        self.pressure = pressure
        self.rotate = 0.0
        self.particles = numberOfParticles
        self.ESCAPE = b'\x1b'  # deneyerek buldum
        self.min_pressure = min(self.pressure)
        self.max_pressure = max(self.pressure)

    def ensure_dir(self, f):
        d = os.path.dirname(f)
        if not os.path.exists(d):
            os.makedirs(d)

    def gnuplot(self, aContainer):

        firstFileName = "results/voropy_v.gnu"
        secondFileName = "results/voropy_p.gnu"
        self.ensure_dir(firstFileName)
        self.ensure_dir(secondFileName)
        verticesFile = open(firstFileName, "w")
        pointsFile = open(secondFileName, "w")

        for cell in aContainer:
            pointsFile.write("%s %s %s %s\n" % (cell.id, cell.pos[0], cell.pos[1], cell.pos[2]))
            for x in range(0, len(cell.vertices())):
                verticesFile.write("%s %s %s\n" % (cell.vertices()[x][0], cell.vertices()[x][1], cell.vertices()[x][2]))
                pass
            verticesFile.write("\n\n\n")
            pass


    def graphr(self, aContainer, timeStep):

        fileName = "results/pressureOut{}.csv".format(timeStep)
        self.ensure_dir(fileName)
        pressureOutFile = open(fileName, "w")

        pressureOutFile.write("DataFormat, 102\n")
        pressureOutFile.write("Memo1\n")
        pressureOutFile.write("Memo2\n")

        for x in range(0, self.particles):
            pressureOutFile.write("%s %s %s\n" %
                                  (aContainer[x].pos[0], aContainer[x].pos[1], self.pressure[x]))
            pass

        pressureOutFile.close()


    def pressureAtParticle(self, particle, timeStep, pressure):

        fileName = "results/pressureAtParticle{}.dat".format(particle)
        self.ensure_dir(fileName)

        pressureAtParticleFile = open(fileName, "a")

        pressureAtParticleFile.write("%s %s\n" % (timeStep, pressure[particle]))

        pressureAtParticleFile.close()


    def pressureAtParticleDetail(self, particle, timeStep, iteration, pressure):

        fileName = "results/pressureAtParticleDetail{}.dat".format(particle)
        self.ensure_dir(fileName)

        pressureAtParticleDetailFile = open(fileName, "a")

        pressureAtParticleDetailFile.write("%s %s %s\n" % (timeStep, iteration, pressure[particle]))

        pressureAtParticleDetailFile.close()


    def permeabilityGraphr(self, aContainer, permeabilityX, permeabilityY, permeabilityZ):

        firstFileName = "results/permXGraphr.csv"
        self.ensure_dir(firstFileName)

        permXGraphrFile = open(firstFileName, "w")

        permXGraphrFile.write("DataFormat, 102\n")
        permXGraphrFile.write("Memo1\n")
        permXGraphrFile.write("Memo2\n")

        for x in range(0, len(permeabilityX)):
            permXGraphrFile.write("%s %s %s\n" %
                                  (aContainer[x].pos[0], aContainer[x].pos[1], permeabilityX[x]))
            pass

        permXGraphrFile.close()

        secondFileName = "results/permYGraphr.csv"
        self.ensure_dir(secondFileName)

        permYGraphrFile = open(secondFileName, "w")

        permYGraphrFile.write("DataFormat, 102\n")
        permYGraphrFile.write("Memo1\n")
        permYGraphrFile.write("Memo2\n")

        for x in range(0, len(permeabilityY)):
            permYGraphrFile.write("%s %s %s\n" %
                                  (aContainer[x].pos[0], aContainer[x].pos[1], permeabilityY[x]))
            pass

        permYGraphrFile.close()

        thirdFileName = "results/permZGraphr.csv"
        self.ensure_dir(thirdFileName)

        permZGraphrFile = open(thirdFileName, "w")

        permZGraphrFile.write("DataFormat, 102\n")
        permZGraphrFile.write("Memo1\n")
        permZGraphrFile.write("Memo2\n")

        for x in range(0, len(permeabilityZ)):
            permZGraphrFile.write("%s %s %s\n" %
                                  (aContainer[x].pos[0], aContainer[x].pos[1], permeabilityZ[x]))
            pass

        permZGraphrFile.close()


    def writeBottomHolePressure(self, particle, timeStep, pressure):

        fileName = "results/pressureWfAtParticle{}.dat".format(particle)
        self.ensure_dir(fileName)

        bottomHolePressureFile = open(fileName, "a")

        bottomHolePressureFile.write("%s %s\n" % (timeStep, pressure))

        bottomHolePressureFile.close()


    def writeDataFrame(self, timestep, data):
        #to write any data frame such as coefficient

        fileName = "results/dataFrame{}.xls".format(timestep)
        self.ensure_dir(fileName)

        dataFrame = pd.DataFrame(data)
        dataFrameFile = open(fileName, "w")

        dataFrameFile.write(dataFrame.to_string())

        dataFrameFile.close()

