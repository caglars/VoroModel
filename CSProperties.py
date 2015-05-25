__author__ = 'caglars'

import CSData
import math

class CSFluidProperties():
    def __init__(self):
        self.reader = CSData.CSDataReader()
        self.referenceDepth = self.reader.readSingleFloatValue("REFDEPTH")
        self.referencePressure = self.reader.readSingleFloatValue("REFPRES")
        self.referenceDensity = self.reader.readSingleFloatValue("REFDENSITY")
        self.liquidCompressibility = self.reader.readSingleFloatValue("LIQUIDCOMPRESSIBILITY")
        self.referenceFVF = self.reader.readSingleFloatValue("REFFVF")
        self.pressureAtBubblePoint = self.reader.readSingleFloatValue("BPPRESSURE")
        self.viscosityAtBubblePoint = self.reader.readSingleFloatValue("BPVISCOSITY")
        #self.referenceDepth = 7000.
        #self.referencePressure = 3031.
        #self.fluidDensity = 62.4

    def findGammaFluid(self, pressure):
        # TODO the value should change with pressure and temperature
        gammaConstant = 0.21584e-3
        gravityAcceleration = 32.17
        # Eqn 2.101 or Eqn 8.88 in Ertekin
        self.fluidDensity = self.referenceDensity*(1+self.liquidCompressibility*(pressure - self.referencePressure))
        print(self.fluidDensity)
        gammaFluid = gammaConstant*self.fluidDensity*gravityAcceleration
        return gammaFluid

    def findFormationVolumeFactor(self, pressure):
        newFormationVolumeFactor = self.referenceFVF/(1+self.liquidCompressibility*(pressure-self.referencePressure))
        return newFormationVolumeFactor

    def findViscosity(self, pressure, correlation):
        # 1 = Vasquez and Beggs Correlation for Undersaturated oil
        if correlation == 1:
            constant1 = 2.6
            constant2 = 1.187
            constant3 = -11.513
            constant4 = -8.98*10**(-5)
            m = constant1*pressure**constant2 * math.exp(constant3+constant4*pressure)
            oilViscosity = self.viscosityAtBubblePoint*(pressure/self.pressureAtBubblePoint)**m
            return oilViscosity