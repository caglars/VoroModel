__author__ = 'caglars'

import pandas as pd

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


    def gnuplot(self, aContainer):


        verticesFile = open("voropy_v.gnu", "w")
        pointsFile = open("voropy_p.gnu", "w")

        for cell in aContainer:
            pointsFile.write("%s %s %s %s\n" % (cell.id, cell.pos[0], cell.pos[1], cell.pos[2]))
            for x in range(0, len(cell.vertices())):
                verticesFile.write("%s %s %s\n" % (cell.vertices()[x][0], cell.vertices()[x][1], cell.vertices()[x][2]))
                pass
            verticesFile.write("\n\n\n")
            pass


    def graphr(self, aContainer, timeStep):

        pressureOutFile = open("pressureOut{0}.csv".format(timeStep), "w")

        pressureOutFile.write("DataFormat, 102\n")
        pressureOutFile.write("Memo1\n")
        pressureOutFile.write("Memo2\n")

        for x in range(0, self.particles):
            pressureOutFile.write("%s %s %s\n" %
                                  (aContainer[x].pos[0], aContainer[x].pos[1], self.pressure[x]))
            pass

        pressureOutFile.close()


    def pressureAtParticle(self, particle, timeStep, pressure):
        pressureAtParticleFile = open("pressureAtParticle{}.dat".format(particle), "a")

        pressureAtParticleFile.write("%s %s\n" % (timeStep, pressure[particle]))

        pressureAtParticleFile.close()


    def pressureAtParticleDetail(self, particle, timeStep, iteration, pressure):
        pressureAtParticleDetailFile = open("pressureAtParticleDetail{}.dat".format(particle), "a")

        pressureAtParticleDetailFile.write("%s %s %s\n" % (timeStep, iteration, pressure[particle]))

        pressureAtParticleDetailFile.close()


    def permeabilityGraphr(self, aContainer, permeabilityX, permeabilityY, permeabilityZ):

        permXGraphrFile = open("permXGraphr.csv", "w")

        permXGraphrFile.write("DataFormat, 102\n")
        permXGraphrFile.write("Memo1\n")
        permXGraphrFile.write("Memo2\n")

        for x in range(0, len(permeabilityX)):
            permXGraphrFile.write("%s %s %s\n" %
                                  (aContainer[x].pos[0], aContainer[x].pos[1], permeabilityX[x]))
            pass

        permXGraphrFile.close()

        permYGraphrFile = open("permYGraphr.csv", "w")

        permYGraphrFile.write("DataFormat, 102\n")
        permYGraphrFile.write("Memo1\n")
        permYGraphrFile.write("Memo2\n")

        for x in range(0, len(permeabilityY)):
            permYGraphrFile.write("%s %s %s\n" %
                                  (aContainer[x].pos[0], aContainer[x].pos[1], permeabilityY[x]))
            pass

        permYGraphrFile.close()

        permZGraphrFile = open("permZGraphr.csv", "w")

        permZGraphrFile.write("DataFormat, 102\n")
        permZGraphrFile.write("Memo1\n")
        permZGraphrFile.write("Memo2\n")

        for x in range(0, len(permeabilityZ)):
            permZGraphrFile.write("%s %s %s\n" %
                                  (aContainer[x].pos[0], aContainer[x].pos[1], permeabilityZ[x]))
            pass

        permZGraphrFile.close()


    def writeBottomHolePressure(self, particle, timeStep, pressure):
        bottomHolePressureFile = open("pressureWfAtParticle{}.dat".format(particle), "a")

        bottomHolePressureFile.write("%s %s\n" % (timeStep, pressure))

        bottomHolePressureFile.close()


    def writeDataFrame(self, timestep, data):
        #to write any data frame such as coefficient

        dataFrame = pd.DataFrame(data)
        dataFrameFile = open("dataFrame{}.xls".format(timestep), "w")

        dataFrameFile.write(dataFrame.to_string())

        dataFrameFile.close()

