library(TAPcodeV2)

# Read in the TDMS or simulation
TAPexperiment = readTAP("../../../Downloads/11011216.tdms")
plot(TAPexperiment$AMU40$time, TAPexperiment$AMU40$pulses[,5])

# apply baseline correction
TAPexperiment$AMU40$options$baselineStart = 0
TAPexperiment$AMU40$options$baselineEnd = .1

TAPexperiment$AMU40$options$timeStart = .1
TAPexperiment$AMU40$options$timeEnd = 2

TAPexperiment$AMU28$options$baselineStart = 0
TAPexperiment$AMU28$options$baselineEnd = .1

TAPexperiment$AMU28$options$timeStart = .1
TAPexperiment$AMU28$options$timeEnd = 2

TAPexperiment$AMU4$options$baselineStart = 3.75
TAPexperiment$AMU4$options$baselineEnd = 4

TAPexperiment$AMU4$options$timeStart = 0
TAPexperiment$AMU4$options$timeEnd = 1.9

TAPexperiment = baselineCorrect(TAPexperiment, c("AMU40", "AMU28", "AMU4") )
plot(TAPexperiment$AMU40$time, TAPexperiment$AMU40$pulses[,5])
lines(TAPexperiment$AMU28$time, TAPexperiment$AMU28$pulses[,5] * 2, col = "red")
# determine moments

TAPexperiment = moments(TAPexperiment, c("AMU40", "AMU28", "AMU4"))
plot(TAPexperiment$AMU40$moments$M0)

# Scale each AMU by a specific calibration coefficient based on previously collected inert zone data
TAPexperiment$AMU28$options$calibrationCoef = 1.25

TAPexperiment$AMU28$pulses = TAPexperiment$AMU28$pulses * TAPexperiment$AMU28$options$calibrationCoef

plot(TAPexperiment$AMU40$time, TAPexperiment$AMU40$pulses[,5])
lines(TAPexperiment$AMU28$time, TAPexperiment$AMU28$pulses[,5], col = "red")

# calculate Calibration coefficient

plot(TAPexperiment$AMU4$moments$M0 / TAPexperiment$AMU40$moments$M0)
plot(TAPexperiment$AMU4$pulses[,5])


# Transform the inert to look like the reactant gas.  (Graham's Law)
names(TAPexperiment)
TAPexperiment = inertGraham(TAPexperiment, "AMU40", 28)
names(TAPexperiment)

# Y-Procedure
TAPexperiment$AMU28$options$inert = "AMU40to28"
TAPexperiment$AMU28$options$isProduct = F

TAPexperiment = gasCon(TAPexperiment, "AMU28", method = "g")
plot(TAPexperiment$AMU28gasCon$pulses[,30])

TAPexperiment = reactionRate(TAPexperiment, "AMU28", method = "g")
plot(TAPexperiment$AMU28rate$pulses[,30])





