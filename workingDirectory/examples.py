### a set of examples
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import TAPas



bedPorosity = 0.5
inertZone1 = .5
inertZone2 = inertZone1
diffusion = 0.5
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
plt.show(block = False)

plt.close()
tapSolverInert = pd.read_csv("test_folder/flux_data/Inert-1.csv",header = None)
tapSolverTime = tapSolverInert.iloc[:,0]
tapSolverInertFlux = tapSolverInert.iloc[:,1]
plt.plot(timeVector, inertFlux)
plt.scatter(tapSolverTime, tapSolverInertFlux, color = "red")
plt.show(block = False)



### Irreversible Rate
plt.close()
tempFlux = TAPas.genIrreversible(30, inertFlux)
plt.plot(timeVector, inertFlux)
plt.plot(timeVector, tempFlux, color = 'red')
plt.show(block=False)

plt.close()
tapSolverReactant = pd.read_csv("test_folder/flux_data/CO.csv",header = None)
tapSolverReactantFlux = tapSolverReactant.iloc[:,1]
plt.plot(timeVector, tempFlux)
plt.scatter(tapSolverTime, tapSolverReactantFlux, color = "red")
plt.show(block=False)

plt.close()
tapSolverReactant = pd.read_csv("test_folder/flux_data/CO.csv",header = None)
tapSolverReactantFlux = tapSolverReactant.iloc[:,1]
plt.plot(timeVector, inertFlux)
plt.plot(timeVector, tempFlux, color = 'green')
plt.scatter(tapSolverTime, tapSolverReactantFlux, color = "red")
plt.show(block=False)




### NO2


plt.close()
tapSolverInert = pd.read_csv("test_folder/flux_data/Inert-1.csv",header = None)
tapSolverTime = tapSolverInert.iloc[:,0]
tapSolverInertFlux = tapSolverInert.iloc[:,1]
tapSolverReactant = pd.read_csv("test_folder/flux_data/AB.csv",header = None)
tapSolverAB = tapSolverReactant.iloc[:,1]
tapSolverReactant = pd.read_csv("test_folder/flux_data/A.csv",header = None)
tapSolverA = tapSolverReactant.iloc[:,1]
tapSolverReactant = pd.read_csv("test_folder/flux_data/B.csv",header = None)
tapSolverB = tapSolverReactant.iloc[:,1]
plt.scatter(tapSolverTime, tapSolverInertFlux - tapSolverAB, color = "red")
plt.scatter(tapSolverTime, tapSolverA, color = "black")
plt.scatter(tapSolverTime, tapSolverB, color = "green")
plt.show(block=False)

plt.close()
tapSolverInert = pd.read_csv("test_folder/flux_data/Inert-1.csv",header = None)
tapSolverTime = tapSolverInert.iloc[:,0]
tapSolverInertFlux = tapSolverInert.iloc[:,1]
tapSolverReactant = pd.read_csv("test_folder/flux_data/AB.csv",header = None)
tapSolverAB = tapSolverReactant.iloc[:,1]
tapSolverReactant = pd.read_csv("test_folder/flux_data/A.csv",header = None)
tapSolverA = tapSolverReactant.iloc[:,1]
tapSolverReactant = pd.read_csv("test_folder/flux_data/B.csv",header = None)
tapSolverB = tapSolverReactant.iloc[:,1]
plt.scatter(tapSolverTime, tapSolverInertFlux)
plt.scatter(tapSolverTime, tapSolverAB, color = "red")
plt.show(block=False)

np.trapz(tapSolverInertFlux, tapSolverTime) - (np.trapz(tapSolverAB, tapSolverTime) + np.trapz(tapSolverA, tapSolverTime) + np.trapz(tapSolverB, tapSolverTime))
np.trapz(tapSolverInertFlux, tapSolverTime) - np.trapz(tapSolverAB, tapSolverTime) - np.trapz(tapSolverA, tapSolverTime)

np.trapz(tapSolverA, tapSolverTime) / np.trapz(tapSolverAB, tapSolverTime)


plt.close()
tapSolverReactant = pd.read_csv("test_folder/flux_data/CO.csv",header = None)
tapSolverReactantFlux = tapSolverReactant.iloc[:,1]
plt.plot(timeVector, inertFlux)
plt.plot(timeVector, tempFlux, color = 'green')
plt.scatter(tapSolverTime, tapSolverReactantFlux, color = "red")
plt.show(block=False)



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