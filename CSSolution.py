__author__ = 'caglars'

import math
import numpy as np


class CSSolver():
    def simpleSolver(self, arraySize, coefficientArray, rightHandSideArray, initialValueArray):
        valueArray = [0 for x in range(arraySize)]
        oldValueArray = [0 for x in range(arraySize)]
        leftHandSideArray = [0 for x in range(arraySize)]

        for i in range(0, arraySize):
            valueArray[i] = initialValueArray[i]
            oldValueArray[i] = valueArray[i]
            pass

        counter = 0

        while True:
            for i in range(0, arraySize):
                LHSIndex = 0
                sumError = 0

                leftHandSideArray[i] = 0

                for j in range(0, arraySize - 1):
                    LHSIndex = (i + j + 1) % arraySize
                    leftHandSideValue = leftHandSideArray[i] \
                                        + coefficientArray[i][LHSIndex] * valueArray[LHSIndex]
                    leftHandSideArray[i] = leftHandSideValue
                    pass
                solution = (rightHandSideArray[i] - leftHandSideArray[i]) / coefficientArray[i][i]
                valueArray[i] = solution
                sumError = sumError + math.fabs(oldValueArray[i] - valueArray[i])
                oldValueArray[i] = valueArray[i]
                pass
            counter = counter + 1
            print("Counter is %s and error is %s" % (counter, sumError))
            if (sumError < 0.0001 or counter > 50):
                break

        return valueArray

    def numpySolver(self, coefficientArray, rightHandSideArray):
        a = np.array(coefficientArray)
        b = np.array(rightHandSideArray)

        return np.linalg.solve(a, b)
