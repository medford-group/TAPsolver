import numpy as np
import pandas as pd
import matplotlib as m_plt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

fig3, ax3 = plt.subplots()

df_p1 = pd.read_csv('./variation_data_time.csv')

x = df_p1.iloc[:,0]

legend_list = ['I','II','III','IV','V','VI']

for k in range(0,len(df_p1.columns)):
	ax3.plot(x,input_list[z],c=color,linestyle=line_in)
