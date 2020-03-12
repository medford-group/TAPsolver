import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import imageio

#def imageGif(pulseNum):

fig,ax = plt.subplots()

CO1 = pd.read_csv("./pulse_3/CO.csv",header=None)
CO2 = pd.read_csv("./pulse_3/*.csv",header=None)
CO3 = pd.read_csv("./pulse_3/O*.csv",header=None)
CO4 = pd.read_csv("./pulse_3/CO*.csv",header=None)
CO5 = pd.read_csv("./pulse_3/r_CO.csv",header=None)

CO6 = pd.read_csv("./pulse_3/O2.csv",header=None)
CO7 = pd.read_csv("./pulse_3/CO2.csv",header=None)


x_axis = []

for k in range(0,150):
	x_axis.append(k*0.001)

#CO1.iloc[:,1] += 3


#max1 = max(CO1.iloc[:150,1])
#max2 = max(CO2.iloc[:150,1])
#max3 = max(CO3.iloc[:150,1])
#max4 = max(CO4.iloc[:150,1])
#max5 = max(CO5.iloc[:150,1])
#max6 = max(CO6.iloc[:150,1])
#max7 = max(CO7.iloc[:150,1])

#for k in range(0,2):
#ax.plot(x_axis,CO1.iloc[:150,1]/max1,label='CO')#
#ax.plot(x_axis,CO2.iloc[:150,1]/max2,label='*')
#ax.plot(x_axis,CO3.iloc[:150,1]/max3,label='O*')
#ax.plot(x_axis,CO4.iloc[:150,1]/max4,label='CO*')
#ax.plot(x_axis,CO5.iloc[:150,1]/max5,label='r_CO')
#ax.plot(x_axis,CO6.iloc[:150,1]/max6,label='O2')
#ax.plot(x_axis,CO7.iloc[:150,1]/max7,label='CO2')

#ax.set_xscale('log')

ax.legend()
fig.canvas.draw()
	
image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	
#return image

plt.show()

#kwargs_write = {'fps':4.0, 'quantizer':'nq'}
#components = list(range(0,2))
#for zoo in range(0,10):
#	components.append(all_steps-1) 

#print(all_steps)
#imageio.mimsave('./output.gif', [imageGif(i) for i in range(1,31)], fps=4)