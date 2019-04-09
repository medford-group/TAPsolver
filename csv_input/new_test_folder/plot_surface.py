import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import csv
import sys
import os
import ufl

pulses = 1

df1 = pd.read_csv('./*.csv',header=None)
df2 = pd.read_csv('./CO*.csv',header=None)
df3 = pd.read_csv('./O2*.csv',header=None)
df4 = pd.read_csv('./O*.csv',header=None)

fig2,ax2 = plt.subplots()
plt.axhline(y=0,color='pink')
ax2.set_xlabel('time (s)')
ax2.set_ylabel('Surface coverage in catalyst zone')

for k in range(1,pulses+1):
	#ax2.plot(df1[0],df1[k],color='k')
	#ax2.plot(df2[0],df2[k],color='b')
	#ax2.plot(df3[0],df3[k],color='r')
	ax2.plot(df4[0],df4[k],color='g')

df1 = pd.read_csv('./O2.csv',header=None)
df2 = pd.read_csv('./CO.csv',header=None)
df3 = pd.read_csv('./CO2.csv',header=None)
df4 = pd.read_csv('./Inert.csv',header=None)

fig3,ax3 = plt.subplots()
plt.axhline(y=0,color='pink')
for k in range(1,pulses+1):
	ax3.plot(df1[0],df1[k],color='k')
	#ax3.plot(df2[0],df2[k],color='b')
	#ax3.plot(df3[0],df3[k],color='r')
	#ax3.plot(df4[0],df4[k],color='g')
plt.show()

plt.show()




#	fig2, ax2 = plt.subplots()
#	ax2.set_xlabel('$t (s)$')
#	if reactor == 'tap':
#		ax2.set_ylabel('$Flux (molecules/s)$')
#	elif reactor == 't_pfr' or 't_pfr_diff':
#		ax2.set_ylabel('$Outlet Concentration (molecules/cm3)$')
#	#zef = np.linspace(dx_r, 1.-dx_r, grid_points+1)
#	
#	legend_label = []
#	header = "time"
#	for k in range(0,gas_phase):
#		legend_label.append(reacs[k])
#		header = header+","+reacs[k]
#	legend_label.append("Inert")
#	header = header+",Inert"
#	colors = ['b','r','g','m','k','y','c']