### a set of examples
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import TAPas



bedPorosity = 0.4
inertZone1 = .018
inertZone2 = inertZone1
diffusion = 0.002
timeStep = .001
lenFlux = 3000
smoothing = 3.0
gasMass = 40.0
inertMass = 40.0
catalystZone = 0

reactorLength = inertZone1 + inertZone2 + catalystZone
timeVector = np.linspace(timeStep, timeStep * lenFlux, num = lenFlux)

### Standard Diffusion Curve
# the base function assumes diffusion coefficient of 0.4, reactor length of 1, timeStep of .001, lenFlux of 3000 (3 seconds), and a bed porosity of 0.4
inertFlux = TAPas.standardDiffusionCurve(diffusion=diffusion, reactorLength = reactorLength, timeStep = timeStep, lenFlux = lenFlux, bedPorosity = bedPorosity)
plt.plot(timeVector, inertFlux)
plt.show()

### Irreversible Rate
tempFlux = TAPas.genIrreversible(5, inertFlux)
plt.plot(timeVector, inertFlux)
plt.plot(timeVector, tempFlux, color = 'red')
plt.show()


### Y-Procedure
# dataPath = "~/Desktop/yProcedureTestData.csv"
# exampleData = pd.read_csv(dataPath)
# gasPulse = exampleData["gasPulse"]
# inertPulse = exampleData['inert']
# isProduct = False


# test = TAPas.yProcedure(gasPulse = gasPulse, gasMass = gasMass, isProduct = isProduct,inertPulse = inertPulse, inertMass = inertMass, timeStep = timeStep, inertZone1 = inertZone1, catalystZone = catalystZone, diffusion = 0.002)
# plt.plot(test['gasConcentration'])
# plt.show()
# plt.plot(test['reactionRate'])
# plt.show()