__author__ = 'caglars'

import CSData

class CSFluidProperties():
    def __init__(self):
        self.reader = CSData.CSDataReader()
        self.referenceDepth = self.reader.readSingleFloatValue("REFDEPTH")
        self.referencePressure = self.reader.readSingleFloatValue("REFPRES")
        self.referenceDensity = self.reader.readSingleFloatValue("REFDENSITY")
        self.liquidCompressibility = self.reader.readSingleFloatValue("LIQUIDCOMPRESSIBILITY")
        #self.referenceDepth = 7000.
        #self.referencePressure = 3031.
        #self.fluidDensity = 62.4

    def findGammaFluid(self, pressure):
        # TODO the value should change with pressure and temperature
        gammaConstant = 0.21584e-3
        gravityAcceleration = 32.17
        self.fluidDensity = self.referenceDensity*(1-self.liquidCompressibility*(pressure - self.referencePressure))
        print(self.fluidDensity)
        gammaFluid = gammaConstant*self.fluidDensity*gravityAcceleration
        return gammaFluid