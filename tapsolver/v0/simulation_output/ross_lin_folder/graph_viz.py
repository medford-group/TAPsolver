import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import math as mp
import sys

fig2,ax2 = plt.subplots()
plt.axhline(y=0,color='pink')

df1 = pd.read_csv('./linA.txt',sep=" ",header=None)
df2 = pd.read_csv('./linB.txt',sep=" ",header=None)
df3 = pd.read_csv('./inertA.txt',sep=" ",header=None)
df4 = pd.read_csv('../ross_inert_folder/Inert.csv',header=None)
ax2.plot(df4[0],df4[1],color='g')
ax2.plot(df1[0],df1[1],color='b')
ax2.plot(df2[0],df2[1],color='r')
#ax2.plot(df3[0],df3[1],color='g')

values = [10,50,100,500,1000,5000]

#for k in values:
#	print(k)

df1n = pd.read_csv('../csv_input_sim/lin_2050_folder/A.csv',header=None)
df2n = pd.read_csv('../csv_input_sim/lin_2050_folder/B.csv',header=None)
df3n = pd.read_csv('../csv_input_sim/lin_2050_folder/Inert.csv',header=None)
#	#ax2.scatter(df1n[0][::5],df1n[1][::5],color='b',s=6)
ax2.plot(df1n[0],df1n[1],color='b',linestyle='--')
#	#ax2.scatter(df2n[0][::5],df2n[1][::5],color='r',s=6)
ax2.plot(df2n[0],df2n[1],color='r',linestyle='--')
#	#ax2.scatter(df3n[0][::5],df3n[1][::5],color='g',s=6)
ax2.plot(df3n[0],df3n[1],color='g',linestyle='--')


df1n = pd.read_csv('../csv_input_sim/lin_4100_folder/A.csv',header=None)
df2n = pd.read_csv('../csv_input_sim/lin_4100_folder/B.csv',header=None)
df3n = pd.read_csv('../csv_input_sim/lin_4100_folder/Inert.csv',header=None)
#	#ax2.scatter(df1n[0][::5],df1n[1][::5],color='b',s=6)
ax2.plot(df1n[0],df1n[1],color='b',linestyle='--')
#	#ax2.scatter(df2n[0][::5],df2n[1][::5],color='r',s=6)
ax2.plot(df2n[0],df2n[1],color='r',linestyle='--')
#	#ax2.scatter(df3n[0][::5],df3n[1][::5],color='g',s=6)
ax2.plot(df3n[0],df3n[1],color='g',linestyle='--')

#Text box and line labels
props = dict(boxstyle='round', facecolor='white', alpha=0.4)
ax2.text(0.89, 0.77, '\n'.join(('Reactions','A -> A*','A* -> B')), transform=ax2.transAxes, fontsize=12,verticalalignment='top', bbox=props,ha='center')
plt.text(0.1, 0.65, 'A',fontsize=14)
plt.text(0.8, 0.7, 'B',fontsize=14)
plt.text(0.17, 2.65, 'Inert',fontsize=14)

##Legend
styles = ['-', '--']
lines = [Line2D([0], [0], color='k', linewidth=3, linestyle=c) for c in styles]
labels = ['Analytical', 'FEniCS']
plt.legend(lines, labels, facecolor='white',title="Method")

ax2.set_xlim(0,1)
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Dimensionless Flux (1/s)',fontsize=12)

plt.show()