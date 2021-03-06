__author__ = 'caglars'

import linecache
import numpy

class CSDataReader():
    def __init__(self):
        self.myDataFile = "dataFile.dat"
        self.particles = 0

    def readDataFor(self, theString):
        lookup = theString
        with open(self.myDataFile) as myFile:
            for num, line in enumerate(myFile, 0):
                # I have added split() to find only the whole word
                if lookup in line.split():
                    #print("Found at line %s" % num)
                    return num
        return 0

    def readParticles(self):
        particleList = []
        start = self.readDataFor("PARTICLES")
        end = self.readDataFor("ENDPARTICLES")
        #myFile = open(self.myDataFile, 'rb', 0)
        #lineCounter = start
        for lineCounter in range(start+2, end+1):
            line = linecache.getline(self.myDataFile, lineCounter)
            #print("line %s start %s end %s lineCounter %s" % (line, start, end, lineCounter))
            lineParticles = line.split()
            particleList.append([lineParticles[0], lineParticles[1], lineParticles[2]])
            linecache.clearcache()

        self.particles = len(particleList)
        return particleList

    def readValues(self, propertyString):
        valueList = []
        start = self.readDataFor(propertyString)
        end = self.readDataFor("END"+propertyString)
        #print ("start %s end %s" % (start, end))

        for lineCounter in range(start+2, end+1):
            line = linecache.getline(self.myDataFile, lineCounter)
            #print("line %s start %s end %s lineCounter %s" % (line, start, end, lineCounter))
            words = line.split()
            #print("words %s" % words)

            for text in words:
                if "*" in text:
                    value = text.split('*')
                    for count in range(0, int(value[0])):
                        valueList.append(float(value[1]))
                    #print("valueList %s" % valueList)
                else:
                    valueList.append(float(text))
            linecache.clearcache()

        return valueList

    def readSingleFloatValue(self, propertyString):
        start = self.readDataFor(propertyString)
        #end = self.readDataFor("END"+propertyString)
        lineCounter = start+2
        line = linecache.getline(self.myDataFile, lineCounter)
        value = float(line)
        return value

    def readSingleIntValue(self, propertyString):
        start = self.readDataFor(propertyString)
        #end = self.readDataFor("END"+propertyString)
        lineCounter = start+2
        line = linecache.getline(self.myDataFile, lineCounter)
        value = int(line)
        return value

    def readWellData(self, dataType):
        #print("readWellRates")
        start = self.readDataFor("WELLS")
        end = self.readDataFor("ENDWELLS")
        #print ("start %s end %s" % (start, end))
        wellData = numpy.zeros(self.particles)
        for lineCounter in range(start+2, end+1):
            #print("here")
            line = linecache.getline(self.myDataFile, lineCounter)
            #print("line %s start %s end %s lineCounter %s" % (line, start, end, lineCounter))
            propertyList = line.split()
            if dataType=="FLOWRATE":
                if propertyList[1] == "FLOWRATE":
                    wellData[int(propertyList[0])] = propertyList[2]
            if dataType=="SKIN":
                if propertyList[1] == "SKIN":
                    wellData[int(propertyList[0])] = propertyList[2]
            if dataType=="RW":
                if propertyList[1] == "RW":
                    wellData[int(propertyList[0])] = propertyList[2]
            if dataType=="PERFTHICK":
                if propertyList[1] == "PERFTHICK":
                    wellData[int(propertyList[0])] = propertyList[2]
            linecache.clearcache()
        return wellData