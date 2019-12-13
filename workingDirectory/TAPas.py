## TAPas Python Functions
import numpy as np

def gProcedure(gasPulse, gasMass, isProduct, inertPulse, inertMass, timeStep, inertZone1, catalystZone, inertZone2=None, diffusion = None, bedPorosity=0.4):
    #' G-Procedure
    #'
    #' Reconstruct concentrations and reaction rates for reactants and products using the G-procedure.
    #' This function uses the Gamma Distribution to translate exit flow into rate and
    #' concentration data from a catalytic process. 
    #'
    #' References: TBD
    #'
    #' Dependencies: Numpy
    #'
    #' Author: M. Ross Kunz
    #'
    ##################
    ################## Definition of inputs
    ##################
    # 
    # gasPulse(numpy array vector): the gas flux
    # gasMass(float): The AMU mass of the gasPulse
    # isProduct(bool): True if gasPulse is a product
    # inertPulse(numpy array vector): the inert flux
    # inertMass(float): The AMU mass of the inertPulse
    # timeStep(float): Time between collection of flux information
    # inertZone1(float): The length of the first inert zone
    # catalystZone(float): The length of the catalyst zone
    # inertZone2(float or None): The length of the second inert zone. If set to None, set to inertZone1
    # diffusion(float or None): The diffusion coefficient of the inert gas. If set to None, will calculate it via moments
    # bedPorosity(float): The bed porosity within the reactor (assuming the same particle size in the inert and catalyst zone)
    # 
    ##################
    ################## Output as dictionary
    ##################
    #
    # gasConcentration(numpy array vector): The calculated gas concentration
    # reactionRate(numpy array vector): The calculated reaction rate
    #
    ##################
    ##################
    ##################
    #
    #
    # set inertZone2 if symmetric reactor
    if inertZone2 is None:
        inertZone2 = inertZone1
    lengthOfInertZone = inertZone1 + inertZone2
    reactorLength = lengthOfInertZone + catalystZone
    # Calculate diffusion if None
    if diffusion is None:
        diffusion = bedPorosity * (inertZone1 + catalystZone + inertZone2)**2 * np.trapz(inertPulse) / (2 * np.trapz(inertPulse * np.arange(0, (len(inertPulse) * timeStep), step=timeStep)))
    lenFlux = len(inertPulse)
    # Calculating Time
    timeVector = np.reshape(np.linspace(timeStep, timeStep * lenFlux, num = lenFlux), (1, lenFlux))
    # Mass time correction if there is a difference between the gas and inert flux mass
    if inertMass != gasMass:
        inertPulse = grahamScaling(inertPulse, inertMass, gasMass, timeStep)
    # Required Moments
    fluxArea = np.trapz(gasPulse, timeVector)
    if isProduct:
        tempReactant = gasPulse
    else: 
        tempReactant =  inertPulse - gasPulse
    # Catalyst Ratio
    catalystRatio = inertZone2 / reactorLength
    ##
    # Gas Concentration
    ##
    # altering the flux wrt the alpha parameter in the gamma distribution
    timeScalar = (1 - catalystRatio**2) * (1/6)
    concentration = gasPulse * timeVector**(-timeScalar) 
    concentration[1] = 0
    # Calculate areas
    concentrationArea = np.trapz(concentration, timeVector) 
    # Appropriately scale the M0
    concentration = concentration / concentrationArea * fluxArea * inertZone2 / diffusion 
    ##
    # Reaction Rate
    ##
    # have to be aware of inert gas being less than reactant gas in initial timing
    # Arbitrarily setting the first 10 values to be zero to account for noise while waiting for initial collection time
    tempReactant[:10] = 0
    # altering the flux wrt the alpha parameter in the gamma distribution
    timeScalar = (1 - catalystRatio**2) * (3 / 2) 
    rate = tempReactant * timeVector**(-timeScalar)
    rate[:10] = 0
    # Appropriately scale the M0
    areaScalar = np.trapz(tempReactant, timeVector) / np.trapz(rate, timeVector)
    rate = rate * areaScalar * timeScalar
    # return result as a dictionary
    result = {
        "gasConcentration":concentration,
        "reactionRate":rate
    }
    return(result)



def genIrreversible(rate, inertPulse):
    #' Inert transformation to irreversible flux
    #'
    #' Scale the standard diffusion curve by the exponential of time and the reaction rate
    #'
    #' Reference: 'TAP-2: An interrogative kinetics approach' by Gleaves et al
    #'
    #' Dependencies: Numpy
    #'
    #' Author: M. Ross Kunz
    #'
    ##################
    ################## Definition of inputs
    ##################
    #
    # rate(float): The rate of the irreversible reaction
    # inertPulse(numpy array vector): the inert flux
    # 
    ##################
    ################## Output as numpy array
    ##################
    #
    # transformedFlux(numpy array vector): Inert flux transformed by the irreversible rate
    #
    ##################
    ##################
    ##################
    #
    #
    # Calculate Graham's Constant
    timeVector = np.linspace(.001, .001 * len(inertPulse), num = len(inertPulse))
    transformedFlux = np.exp(-rate * timeVector) * inertPulse
    return(transformedFlux)


def grahamScaling(inertPulse, inertMass, gasMass, timeStep):
    #' Flux transformation by Graham's Law
    #'
    #' Create a new flux based on the time manipulation of Graham's Law
    #'
    #' Reference: TBD M. Ross Kunz
    #'
    #' Dependencies: Numpy
    #'
    #' Author: M. Ross Kunz
    #'
    ##################
    ################## Definition of inputs
    ##################
    #
    # inertPulse(numpy array vector): the inert flux
    # inertMass(float): The AMU mass of the inertPulse
    # gasMass(float): The AMU mass of the non-reactant pulse
    # timeStep(float): Time between collection of flux information
    # 
    ##################
    ################## Output as numpy array
    ##################
    #
    # transformedFlux(numpy array vector): Inert flux transformed by mass to determine diffusion response of new gas
    #
    ##################
    ##################
    ##################
    #
    #
    # Calculate Graham's Constant
    grahamsConstant = np.sqrt(inertMass / gasMass)
    roundPt = int(np.floor(len(inertPulse) / grahamsConstant))
    # Calculate Time
    timeVector = np.linspace(0, timeStep * (len(inertPulse) - 1),num = len(inertPulse))
    # Transform the time and flux
    if grahamsConstant >= 1:
        transformedFlux = np.interp(x = np.linspace(0, max(timeVector), num = roundPt), xp = timeVector, fp = inertPulse)
        trailingZeros = np.repeat(transformedFlux[-1], len(timeVector) - roundPt)
        transformedFlux = np.concatenate((transformedFlux, trailingZeros))
    else:
        transformedFlux = np.interp(x = timeVector * grahamsConstant, xp = timeVector, fp = inertPulse)
    # Normalizing back to the original M0
    transformedFlux = transformedFlux / np.trapz(transformedFlux) * np.trapz(inertPulse) 
    return(transformedFlux)


def reactivities(reactantPulse, reactantMass, inertPulse, inertMass, timeStep, inertZone1, catalystZone, productPulse = None, productMass = None, inertZone2 = None, bedPorosity = 0.4):
    #' Reactivities
    #'
    #' Create moments and reactivities for inert, reactant and products
    #'
    #' Reference: 'Precise non-steady state characterization of solid active materials with no preliminary mechanistic assumptions' by Constales et al.
    #'
    #' Dependencies: Numpy
    #'
    #' Author: M. Ross Kunz
    #'
    ##################
    ################## Definition of inputs
    ##################
    #
    # reactantPulse(numpy array vector): the gas flux
    # reactantMass(float): The AMU mass of the reactantPulse
    # inertPulse(numpy array vector): the inert flux
    # inertMass(float): The AMU mass of the inertPulse
    # timeStep(float): Time between collection of flux information
    # inertZone1(float): The length of the first inert zone
    # catalystZone(float): The length of the catalyst zone
    # productPulse(numpy array vector): the product flux
    # productMass(float): The AMU mass of the productPulse
    # inertZone2(float or None): The length of the second inert zone. If set to None, set to inertZone1
    # bedPorosity(float): The bed porosity within the reactor (assuming the same particle size in the inert and catalyst zone)
    # 
    ##################
    ################## Output as dictionary of moments and reactivites
    ##################
    #
    # inert(dictionary): M0, M1, diffusion and residenceTime of the inert
    # reactant(dictionary): M0, M1, M2, diffusion, residenceTime, R0, R1, R2 of the reactant
    # product(dictionary): M0, M1, M2, diffusion, residenceTime, R0, R1, R2 of the product
    #
    ##################
    ##################
    ##################
    #
    #
    # Reactor parameters
    if inertZone2 is None:
        inertZone2 = inertZone1
    lengthOfInertZone = inertZone1 + inertZone2
    reactorLength = lengthOfInertZone + catalystZone
    lenFlux = len(inertPulse)
    # Calculating Time
    timeVector = np.reshape(np.linspace(timeStep, timeStep * lenFlux, num = lenFlux), (1, lenFlux))
    # Required Moments
    inertM0 = np.trapz(inertPulse, timeVector)
    inertM1 = np.trapz(inertPulse * timeVector, timeVector)
    reactantM0 = np.trapz(reactantPulse, timeVector) / inertM0
    reactantM1 = np.trapz(reactantPulse * timeVector, timeVector) / inertM0
    reactantM2 = np.trapz(reactantPulse * timeVector**2, timeVector) / inertM0
    # Diffusion Coefficients
    inertDiffusion = bedPorosity * reactorLength**2 * inertM0 / (2 * inertM1)
    reactantDiffusion = inertDiffusion * np.sqrt(inertMass / reactantMass)
    # Residence Times
    residenceTimeInertZone = catalystZone * lengthOfInertZone / (2 * reactantDiffusion)
    residenceTimeCatZone = bedPorosity * lengthOfInertZone**2 / (2 * reactantDiffusion)
    # Dictionary of Inert Results
    inert = {
        "M0":inertM0,
        "M1":inertM1,
        "diffusion":inertDiffusion,
        "residenceTime":residenceTimeInertZone
    }
    # Convenient pre calculations
    residenceTimeRatio = residenceTimeInertZone / residenceTimeCatZone
    reactantM1ratio = reactantM1 / reactantM0
    residenceTimereactantM0 = 1 / residenceTimeCatZone * reactantM0
    # Reactivities for the Reactant
    reactantR0 = -1 / residenceTimeCatZone + residenceTimereactantM0
    reactantR1 = -2 / 3 * residenceTimeRatio - 1 / 3 * residenceTimeInertZone * residenceTimereactantM0 + residenceTimeRatio * residenceTimereactantM0
    reactantR2 = 4 / 45 * residenceTimeInertZone * residenceTimeRatio + 7 / 90 * residenceTimeInertZone**2 * residenceTimereactantM0 - 1 / 3 * residenceTimeInertZone * reactantM1ratio * residenceTimereactantM0 - 1 / 2 * reactantM2 * residenceTimereactantM0 / reactantM0 + reactantM1ratio**2 * residenceTimereactantM0
    # Dictionary of Reactant Results
    reactant = {
        "M0":reactantM0,
        "M1":reactantM1,
        "M2":reactantM2,
        "diffusion":reactantDiffusion,
        "residenceTime":residenceTimeCatZone,
        "R0":reactantR0,
        "R1":reactantR1,
        "R2":reactantR2
    }
    # Product information
    product = {}
    if productPulse is not None:
        # Product Moments
        productM0 = np.trapz(productPulse, timeVector) / inertM0
        productM1 = np.trapz(productPulse * timeVector, timeVector) / inertM0
        productM2 = np.trapz(productPulse * timeVector**2, timeVector) / inertM0
        # Diffusion and Residence Time
        productDiffusion = inertDiffusion * np.sqrt(inertMass / productMass)
        residenceTimeProduct = catalystZone * lengthOfInertZone / (2 * productDiffusion)
        productR0 = productM0 / (residenceTimeProduct * reactantM0)
        productR1 = productR0 * (residenceTimeCatZone / 12 * (8 * reactantM0 + 3 + 9 * reactantDiffusion / productDiffusion) + residenceTimeProduct * reactantM0 * reactantR1 - productM1 / productM0)
        productR2line1 = productM2 / productM0 - 19 / 16 * reactantDiffusion**2 * residenceTimeInertZone**2 / productDiffusion**2 - 19 / 16 * reactantDiffusion * residenceTimeInertZone / productDiffusion * (productR0 * ((3 + 8 * reactantM0) * residenceTimeInertZone + 12 * reactantM0 * reactantR1 * residenceTimeProduct) - 12 * productR1)
        productR2line2 = 1 / 6 * productR1 / productR0 * ((3 + 8 * reactantM0) * residenceTimeInertZone + 12 * reactantM0 * reactantR1 * residenceTimeProduct)
        productR2line3 = -residenceTimeInertZone**2 * (5 / 48 + 1 / 45 * reactantM0 * (23 + 40 * reactantM0)) 
        productR2line4 = -reactantM0 * (3 + 16 * reactantM0) * reactantR1 * residenceTimeProduct * residenceTimeInertZone - 2 * reactantM0 * residenceTimeProduct * (reactantM0 * reactantR1**2 * residenceTimeProduct - reactantR2)
        productR2 = 1 / 2 * productR0 * (productR2line1 + productR2line2 + productR2line3 + productR2line4)
        # Dictionary of Product Results
        product = {
            "M0":productM0,
            "M1":productM1,
            "M2":productM2,
            "diffusion":productDiffusion,
            "residenceTime":residenceTimeProduct,
            "R0":productR0,
            "R1":productR1,
            "R2":productR2
        }
    result = {
        "inert":inert,
        "reactant":reactant,
        "product":product
    }
    return(result)


def standardDiffusionCurve(diffusion = .4, reactorLength = 1, timeStep = .001, lenFlux = 3000, bedPorosity = 0.4):
    #' Standard Diffusion Curve
    #'
    #' Create the inert outlet flux based on the diffusion coefficient
    #'
    #' Dependencies: Numpy
    #'
    #' Reference: 'TAP-2: An interrogative kinetics approach' by Gleaves et al
    #'
    #' Author: M. Ross Kunz
    #'
    ##################
    ################## Definition of inputs
    ##################
    #
    # diffusion(float or None): The diffusion coefficient of the inert gas
    # timeStep(float): Time between collection of flux information
    # lenFlux(float): Number of measured points within the flux
    # reactorLength(float): Length of the reactor
    # bedPorosity(float): The bed porosity within the reactor (assuming the same particle size in the inert and catalyst zone)
    #
    ##################
    ################## Output as 1d numpy array
    ##################
    #
    # sdc(numpy array vector): The standard diffusion curve
    #
    ##################
    ##################
    ##################
    #
    #
    # residence time of the diffusion curve
    mRatio = bedPorosity * reactorLength**2 / (2 * diffusion)
    # create the time associated with the flux as a vector
    timeVector = np.reshape(np.linspace(timeStep, timeStep * lenFlux, num = lenFlux), (1, lenFlux))
    # The number of terms for the summation.  This is originally done with 2000, but not necessary.
    nTerms = 100 
    jTerms = np.linspace(start = 0, stop = nTerms, num = (nTerms + 1))
    #(-1)^n iteration
    negativeIteration = np.repeat(1,(nTerms + 1))
    negativeIteration[jTerms % 2 != 0] = -1
    negWithJ = negativeIteration * (2 * jTerms + 1)
    jTerms = np.reshape(jTerms, ((nTerms + 1), 1))
    # exponential terms within SDC
    exponentTerms = np.matmul(-0.5 * np.pi**2 * (jTerms + 0.5)**2 / mRatio, timeVector)
    exponentWithTime =  np.matmul(negWithJ, np.exp(exponentTerms))
    result = np.pi * exponentWithTime.flatten() / (mRatio * 2)
    return(result)


def yProcedure(gasPulse, gasMass, isProduct, inertPulse, inertMass, timeStep, inertZone1, catalystZone, inertZone2=None, smoothing=3.0, diffusion=None, bedPorosity=0.4):
    #' Y-Procedure
    #'
    #' Reconstruct concentrations and reaction rates for reactants and products using the Y-procedure.
    #' This function uses the Y-Procedure to translate exit flow into rate and
    #' concentration data from a catalytic process. The function takes as input
    #' exit flow measured in pulses of inert, product, and reactant gas (in AMUs),
    #' physical properties of the reactor, and a smoothing parameter sigma. Sigma is used
    #' for simple exponential dampening of the high frequency noise. Depending on the signal-to-noise ratio,
    #' should be varied 2-4. Oversmoothing can be as detrimental as undersmoothing.
    #'
    #' References: 'The Y-procedure: How to extract the chemical transformation rate from reaction-diffusion data with no assumptions on the kinetic model' by Yablonsky et al
    #' References: 'The Y-procedure methodology for the interpretation of transient kinetic data: Analysis of irreversible adsorption' by Redekop et al
    #'
    #' Dependencies: Numpy
    #'
    #' Author: M. Ross Kunz
    #'
    ##################
    ################## Definition of inputs
    ##################
    # 
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
    # 
    ##################
    ################## Output as dictionary
    ##################
    #
    # gasConcentration(numpy array vector): The calculated gas concentration
    # reactionRate(numpy array vector): The calculated reaction rate
    #
    ##################
    ##################
    ##################
    #
    #
    # set inertZone2 if symmetric reactor
    if inertZone2 is None:
        inertZone2 = inertZone1
    # Calculate diffusion if None
    if diffusion is None:
        diffusion = bedPorosity * (inertZone1 + catalystZone + inertZone2)**2 * np.trapz(inertPulse) / (2 * np.trapz(inertPulse * np.arange(0, (len(inertPulse) * timeStep), step=timeStep)))
    # Mass time correction if there is a difference between the gas and inert flux mass
    if inertMass != gasMass:
        inertPulse = grahamScaling(inertPulse, inertMass, gasMass, timeStep)
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
        "gasConcentration":gasConcentration,
        "reactionRate":reactionRate
    }
    return(result)