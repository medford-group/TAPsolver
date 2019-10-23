def yProcedure(gasPulse, gasMass, isProduct, inertPulse, inertMass, timeStep, inertZone1, catalystZone, inertZone2=None, smoothing=3.0, diffusion=None, bedPorosity=0.4):
    #' Y-Procedure
    #'
    #' Reconstruct concentrations and reaction rates for reactants and products using the Y-procedure
    #'
    #' This function uses the Y-Procedure to translate exit flow into rate and
    #' concentration data from a catalytic process. The function takes as input
    #' exit flow measured in pulses of inert, product, and reactant gas (in AMUs),
    #' physical properties of the reactor, and a smoothing parameter sigma. Sigma is used
    #' for simple exponential dampening of the high frequency noise. Depending on the signal-to-noise ratio,
    #' should be varied 2-4. Oversmoothing can be as detrimental as undersmoothing.
    #'
    #' Depends on Numpy
    #'
    #' Author: M. Ross Kunz
    #'
    ##################
    ################## Definition of inputs
    ##################
    # gasPulse(numpy array vector): the gas flux
    # gasMass(float): The AMU mass of the gasPulse
    # isProduct(bool): True if gasPulse is a product
    # inertPulse(numpy array vector): the inert flux
    # inertMass(float): The AMU mass of the inertPulse
    # timeStep(float): Time between collection of flux information
    # inertZone1(float): The length of the first inert zone
    # catalystZone(float): The length of the catalyst zone
    # inertZone2(float or None): The length of the second inert zone. If set to None, set to inertZone1
    # smoothing(float): The amount of exponential smoothing in the Y-Procedure (not well defined in paper, see paragraph 4 section 5.3)
    # diffusion(float or None): The diffusion coefficient of the inert gas. If set to None, will calculate it via moments
    # bedPorosity(float): The bed porosity within the reactor (assuming the same particle size in the inert and catalyst zone)
    ##################
    ################## Output as dictionary
    ##################
    # gasConcentration(numpy array vector): The calculated gas concentration
    # reactionRate(numpy array vector): The calculated reaction rate
    ##################
    ##################
    ##################
    # set inertZone2 if symmetric reactor
    if inertZone2 is None:
        inertZone2 = inertZone1
    # Calculate diffusion if None
    if diffusion is None:
        diffusion = bedPorosity * (inertZone1 + catalystZone + inertZone2)**2 * np.trapz(inertPulse) / (2 * np.trapz(inertPulse * np.arange(0, (len(inertPulse) * timeStep), step=timeStep)))
    # Mass time correction if there is a difference between the gas and inert flux mass
    if inertMass != gasMass:
        grahamsConstant = np.sqrt(inertMass / gasMass)
        roundPt = int(np.floor(len(inertPulse) / grahamsConstant))
        timeVector = np.linspace(0, timeStep * (len(inertPulse) - 1),num = len(inertPulse))
        if grahamsConstant >= 1:
            newFlux = np.interp(x = np.linspace(0, max(timeVector), num = roundPt), xp = timeVector, fp = inertPulse)
            trailingZeros = np.repeat(newFlux[-1], len(timeVector) - roundPt)
            newFlux = np.concatenate((newFlux, trailingZeros))
        else:
            newFlux = np.interp(x = timeVector * grahamsConstant, xp = timeVector, fp = inertPulse)
        newFlux = newFlux / np.trapz(newFlux) * np.trapz(inertPulse) # normalizing back to the original M0
        inertPulse =  newFlux
    # Initial Values
    n = len(gasPulse)
    tau1 = bedPorosity * inertZone1**2 / diffusion
    tau3 = bedPorosity * inertZone2**2 / diffusion
    gamma1 = diffusion / inertZone1
    gamma3 = diffusion / inertZone2
    # Building Omega
    k1 = np.linspace(0.0, n/2, int(n/2 + 1))
    k2 = np.linspace(-1.0, -(n/2 - 1), int(n / 2 - 1))
    k = np.concatenate((k1, k2))
    # Vector of frequencies for signal measured at equally spaced time points.
    # The signal is the experimental exit flow as measured by QMS.
    # Define variables iwt3, iwt1 for use in cFreqDomain, rFreqDomain
    omega = 2 * np.pi * k / (n * timeStep)
    omega[0] = 1.0e-10
    iwt1 = np.sqrt(1j * omega * tau1)
    iwt3 = np.sqrt(1j * omega * tau3)
    smoothScalar = np.exp(-(omega * timeStep * smoothing)**2 / 2)
    # Calculate area i.e. M0
    gasArea = np.trapz(gasPulse)
    inertArea = np.trapz(inertPulse)
    # Calculate Gas Concentration
    if isProduct:
        gasConcentration = np.zeros(len(gasPulse))
    else:
        cFreqDomain = np.sinh(iwt1)/(gamma1 * iwt1)
        cFreqDomain[0] = 1/gamma1
        cFreqDomain = cFreqDomain * smoothScalar
        gasConcentration = np.real(np.fft.ifft(cFreqDomain * np.fft.fft(gasPulse)))
        gasConcentration = gasConcentration / np.trapz(gasConcentration) * gasArea * inertZone2 / diffusion * 2
    # Calculate the Reaction Rate
    rFreqDomain = (np.cosh(iwt3) + np.sqrt((tau1 * gamma1**2) / (tau3 * gamma3**2)) * np.sinh(iwt1) * np.sinh(iwt3) / np.cosh(iwt1)) * smoothScalar
    if isProduct:
        tempReactant = gasPulse
    else: 
        tempReactant =  inertPulse - gasPulse
    # Make sure the rate has the appropriate area
    reactionRate = np.real(np.fft.ifft(rFreqDomain * np.fft.fft(tempReactant)))
    conversion = np.trapz(tempReactant) / inertArea
    conversion = np.trapz(gasPulse) / np.trapz(reactionRate) * conversion / (1 - conversion)
    reactionRate = reactionRate * conversion * 2
    # return result as a dictionary
    result = {
        "gasConcentration" : gasConcentration,
        "reactionRate" : reactionRate
    }
    return(result)


### Example

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
dataPath = "~/Desktop/yProcedureTestData.csv"
exampleData = pd.read_csv(dataPath)
gasPulse = exampleData["gasPulse"]
inertPulse = exampleData['inert']
isProduct = False
bedPorosity = 0.4
inertZone1 = .018
inertZone2 = None
diffusion = 0.002
timeStep = .001
smoothing = 3.0
gasMass = 40.0
inertMass = 40.0
catalystZone = 0

test = yProcedure(gasPulse = gasPulse, gasMass = gasMass, isProduct = isProduct,inertPulse = inertPulse, inertMass = inertMass, timeStep = timeStep, inertZone1 = inertZone1, catalystZone = catalystZone, diffusion = 0.002)
plt.plot(test['gasConcentration'])
plt.show()
plt.plot(test['reactionRate'])
plt.show()