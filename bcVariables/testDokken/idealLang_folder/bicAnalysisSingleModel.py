import time
import matplotlib.pyplot as plt
import imageio
import pandas as pd
import numpy as np
import math
import os
import sys

listValues = ['1234','123','124','134','234','12','13','14','23','24','23','24','34','1','2','3','4','']
#listValues = ['1234','123','124','134','234','12','13','14','23','24','23','24','34','1','2','3','4','']

x = []
y = []

for kstep in listValues:
	
	print(kstep)
	if kstep == '':
		data_set = 'idealLang_fittingBoth_folder/'
	else:
		data_set = 'idealLang_fittingBoth_'+kstep+'_folder/'
	k = 8-len(kstep)
	##time.sleep(1)
	x.append(k)

	n = 6000

	df1CO = pd.read_csv('./flux_data/CO.csv',header=None)
	df2CO = pd.read_csv('./../'+data_set+'/flux_data/CO.csv',header=None)
	df1O2 = pd.read_csv('./flux_data/O2.csv',header=None)
	df2O2 = pd.read_csv('./../'+data_set+'/flux_data/O2.csv',header=None)
	df1CO2 = pd.read_csv('./flux_data/CO2.csv',header=None)
	df2CO2 = pd.read_csv('./../'+data_set+'/flux_data/CO2.csv',header=None)

	#df3 =   (df1CO2.iloc[:,1] - df2CO2.iloc[:,1])**2/df1CO2.iloc[:,1] + (df1CO.iloc[:,1] - df2CO.iloc[:,1])**2/df1CO.iloc[:,1] + (df1O2.iloc[:,1] - df2O2.iloc[:,1])**2/df1O2.iloc[:,1] 
	df3 =   (df1CO2.iloc[:,1] - df2CO2.iloc[:,1])**2 + (df1CO.iloc[:,1] - df2CO.iloc[:,1])**2 + (df1O2.iloc[:,1] - df2O2.iloc[:,1])**2 

	test = df3.sum()


	#print(n*math.log(test/n) + 2*k)
	#print(k*math.log(n) + n*math.log(test/n))
	print(test/n)
	y.append(test/n)
	#y.append(k*math.log(n) + n*math.log(test/n))

fig,ax = plt.subplots()

ax.scatter(x,y,color = 'r', marker='s')

ax.grid( linestyle='-', linewidth=0.25)

plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

ax.tick_params(axis='x', which='minor', bottom=False)

for i, txt in enumerate(listValues):
    ax.annotate(txt, (x[i], y[i]),fontsize='12')

#plt.grid(True)
plt.xticks([4,5,6,7,8])

ax.set_ylim(2e-4,3e-4)

ax.set_xlabel('Number of Parameters')
ax.set_ylabel('BIC Value')

plt.show()