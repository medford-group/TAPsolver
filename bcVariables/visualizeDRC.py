import mpmath
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp

reactions = ['0','1','2']

gasSpecies = ['CO','O2','CO2']

fig2, ax2 = plt.subplots() 

for j in gasSpecies:
	#fig2, ax2 = plt.subplots() 

	#for k in reactions:
	#	#user_data = pd.read_csv('./test_folder/sensitivity/Ga'+k+'/dc_'+j+'.csv',header=None)
	#	user_data = pd.read_csv('./test_folder/RRM_analysis/Ga'+k+'/thinValue/dr_'+j+'.csv',header=None)

	#	if k == '0':
	#		new = user_data[0]
	#	else:
	#		new += user_data[0]
	#	#plt.plot(user_data[0],label=k)

	user_data = pd.read_csv('./test_folder/thin_data/r_'+j+'.csv',header=None)
	#plt.plot(new,label='Total')
	plt.plot(user_data[0],label=j)

	#ax2.axvline(x=abs(user_data[0]).idxmax(), ymin=-1.25, ymax=1.25)

ax2.set_ylabel('Rate (nmol/s)')
ax2.set_xlabel('time')
ax2.set_title('Rate of Formation/Consumption of Species')
ax2.legend()
#ax2.set_ylim(-1,1)
plt.show()

#fig2, ax2 = plt.subplots() 
#ax2.set_ylabel('KDRC')
#ax2.set_xlabel('time')

#for j in gasSpecies:
#	user_data = pd.read_csv('./test_folder/thin_data/r_'+j+'.csv',header=None)
#	plt.plot(user_data[0]/abs(user_data[0]).max(),label=k)

#plt.show()