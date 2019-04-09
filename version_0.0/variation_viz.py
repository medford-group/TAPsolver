import numpy as np
import pandas as pd
import matplotlib as m_plt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

cmap = plt.get_cmap('rainbow')
legend_label = ['I','II','III','IV','V','VI']
df = pd.read_csv('./variation_data_mesh.csv')

fig3, ax3 = plt.subplots()
#ax3.legend(custom_lines,molecules,title="Gas Species")
ax3.set_xlabel('Mesh Size')
ax3.set_ylabel('t (s)')
x = df.iloc[:,0]

#ax3.legend(custom_lines,molecules,title="Gas Species")

for k in range(1,len(df.columns)):
	color = cmap(float(k)/len(df.columns))
	ax3.plot(x,df.iloc[:,k],c=color,label=legend_label[k-1])
ax3.legend(title="Mechanism")
plt.savefig('./mesh_variation.png')
plt.show()

df = pd.read_csv('./variation_data_mesh.csv')

fig4, ax4 = plt.subplots()
#ax3.legend(custom_lines,molecules,title="Gas Species")
ax4.set_xlabel('Time Steps (s)')
ax4.set_ylabel('t (s)')
x = df.iloc[:,0]

#ax3.legend(custom_lines,molecules,title="Gas Species")

for k in range(1,len(df.columns)):
	color = cmap(float(k)/len(df.columns))
	ax4.plot(x,df.iloc[:,k],c=color,label=legend_label[k-1])
ax4.legend(title="Mechanism")
plt.savefig('./time_variation.png')
plt.show()