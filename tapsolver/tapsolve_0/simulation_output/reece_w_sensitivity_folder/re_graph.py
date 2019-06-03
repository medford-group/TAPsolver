import numpy as np
import pandas as pd
import matplotlib as m_plt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
import sys

df_1 = pd.read_csv('./CO.csv')
df_2 = pd.read_csv('./CO2.csv')
df_3 = pd.read_csv('./Inert.csv')
color = ['b','r','c']
legend_label = ['CO','CO2','Inert']
fig4, ax4 = plt.subplots()
#plt.title(r'{\fontsize{30pt}{3em}\selectfont{}{Mean WRFv3.5 LHF\r}{\fontsize{18pt}{3em}\selectfont{}(September 16 - October 30, 2012)}')
#ax4.set_ylabel(r'$\frac{ F }{ N } $ (1/s)',fontsize=15)
ax4.set_ylabel( 'F/Np (1/s)')
ax4.set_xlabel('t (s)')

ax4.plot(df_1.iloc[:,0],df_1.iloc[:,1],c=color[0],label=legend_label[0])
ax4.plot(df_1.iloc[:,0],df_2.iloc[:,1],c=color[2],label=legend_label[1])
ax4.plot(df_1.iloc[:,0],df_3.iloc[:,1],c=color[1],label=legend_label[2])

for k in range(2,len(df_1.columns)):
	ax4.plot(df_1.iloc[:,0],df_1.iloc[:,k],c=color[0])
	ax4.plot(df_1.iloc[:,0],df_3.iloc[:,k],c=color[1])
	ax4.plot(df_1.iloc[:,0],df_2.iloc[:,k],c=color[2])

blue_patch = mpatches.Patch(color='blue', label='CO')
red_patch = mpatches.Patch(color='red', label='CO2')
cyan_patch = mpatches.Patch(color='cyan', label='Inert')

ax4.legend(title="Gas Species",handles=[blue_patch,red_patch,cyan_patch])

ax4.set_xlim([0,0.4])
plt.savefig('./lang_eley_eluc_2.png')
plt.show()