__author__ = 'caglars'

class CSFluidProperties():
    def __init__(self):
        self.referenceDepth = 7000.
        self.referencePressure = 3031.
        self.fluidDensity = 62.4

    def findGammaFluid(self):
        # TODO the value should change with pressure and temperature
        gammaConstant = 0.21584e-3
        gravityAcceleration = 32.17
        gammaFluid = gammaConstant*self.fluidDensity*gravityAcceleration
        return gammaFluid