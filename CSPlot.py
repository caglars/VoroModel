__author__ = 'caglars'

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

