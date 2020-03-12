import matplotlib.pyplot as plt

################################################################

fig,ax = plt.subplots()

meshSize = [100,200,400,800,1000,1200,1400]
meshSizeTime = [5.688,8.6,14.402,28.766,31.642,37.552,41.99]

#name = 'meshSize'
#ax.plot(meshSize,meshSizeTime)
#ax.scatter(meshSize,meshSizeTime,marker = 's')
#ax.set_xlabel('Mesh Size',fontsize=15)
#ax.set_ylabel('Simulation Time (s)',fontsize=15)
#ax.grid(color='k', linestyle='-', linewidth=0.3, alpha=0.5)

################################################################

timeSteps = [1000,1200,1400,1600,1800,2000]
timeStepsTime = [9.91,11.567,13.584,15.334,17.135,19.099]

#name = 'timeSteps'
#ax.plot(timeSteps,timeStepsTime)
#ax.scatter(timeSteps,timeStepsTime,marker = 's')
#ax.set_xlabel('Number of Time Steps',fontsize=15)
#ax.set_ylabel('Simulation Time (s)',fontsize=15)
#ax.grid(color='k', linestyle='-', linewidth=0.3, alpha=0.5)

################################################################

meshRefine = [0, 1, 2, 3, 4, 5, 6, 7]
meshRefineTime = [15.009,14.731,15.214,16.415,19.187,23.209,32.528,50.408]

#name = 'meshRefine'
#ax.plot(meshRefine,meshRefineTime)
#ax.scatter(meshRefine,meshRefineTime,marker = 's')
#ax.set_xlabel('Number of Mesh Refinements',fontsize=15)
#ax.set_ylabel('Simulation Time (s)',fontsize=15)
#ax.grid(color='k', linestyle='-', linewidth=0.3, alpha=0.5)

################################################################

adjointParams = [1,2,3,4,5,6]

adjointParamsAdj1 = [571,669.96,766.68,875.4,975.24,1053.84]
adjointParamsSim1 = [35.815,35.442,37.181,35.528,35.568,37,37.034]

#name = 'adjoint1'
#ax.plot(adjointParams,adjointParamsAdj1)
#ax.scatter(adjointParams,adjointParamsAdj1,marker = 's')
#ax.set_xlabel('Number of Parameters',fontsize=15)
#ax.set_ylabel('Evaluation Time (s)',fontsize=15)
#ax.grid(color='k', linestyle='-', linewidth=0.3, alpha=0.5)

################################################################

adjointParamsSim2 = [11.924,12.024,11.879,11.996,12.364,11.908,11.994]
adjointParamsAdj2 = [174.3,200.28,230.16,236.46,289.2,325.2]

name = 'adjoint2'
ax.plot(adjointParams,adjointParamsAdj2)
ax.scatter(adjointParams,adjointParamsAdj2,marker = 's')
ax.set_xlabel('Number of Parameters',fontsize=15)
ax.set_ylabel('Evaluation Time (s)',fontsize=15)
ax.grid(color='k', linestyle='-', linewidth=0.3, alpha=0.5)

fig.savefig('./'+name+'.png')
plt.show()