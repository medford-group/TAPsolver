import numpy as np
import pandas as pd
import matplotlib as m_plt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

cmap = plt.get_cmap('rainbow')

custom_lines = [Line2D([0],[0],linestyle='--',color='k'), Line2D([0],[0],linestyle='-',color='k')]

constants = ['Ke0','Kd0','Ke1','Ke2','Ke3']
molecules = ['CO','CO2']
pulse_intensity = 2*1.6727777777777779e+19
folder_location = './'

for k in constants:
	fig3, ax3 = plt.subplots()
	output_file_name = k+'_sens'
	for j in molecules:
		gas = j
		parameter = k
		
		df_p1 = pd.read_csv(folder_location+gas+'/pulse_1.csv')
#		df_p2 = pd.read_csv(folder_location+gas+'/pulse_2.csv')
#		df_p3 = pd.read_csv(folder_location+gas+'/pulse_3.csv')
#		df_p4 = pd.read_csv(folder_location+gas+'/pulse_4.csv')
#		df_p5 = pd.read_csv(folder_location+gas+'/pulse_5.csv')
#		df_p6 = pd.read_csv(folder_location+gas+'/pulse_6.csv')
#		df_p7 = pd.read_csv(folder_location+gas+'/pulse_7.csv')
#		df_p8 = pd.read_csv(folder_location+gas+'/pulse_8.csv')
#		df_p9 = pd.read_csv(folder_location+gas+'/pulse_9.csv')
#		df_p10 = pd.read_csv(folder_location+gas+'/pulse_10.csv')
		
		x = df_p1['# t']
		y1 = df_p1[parameter]/pulse_intensity
#		y2 = df_p2[parameter]/pulse_intensity
#		y3 = df_p3[parameter]/pulse_intensity
#		y4 = df_p4[parameter]/pulse_intensity
#		y5 = df_p5[parameter]/pulse_intensity
#		y6 = df_p6[parameter]/pulse_intensity
#		y7 = df_p7[parameter]/pulse_intensity
#		y8 = df_p8[parameter]/pulse_intensity
#		y9 = df_p9[parameter]/pulse_intensity
#		y10 = df_p10[parameter]/pulse_intensity
		if j == 'CO':
			line_in = '--'
		else:
			line_in = '-'
		input_list = [y1]#,y2,y3,y4,y5,y6,y7,y8,y9,y10

		for z in range(0,1):
			color = cmap(float(z)/10)
			ax3.plot(x,input_list[z],c=color,linestyle=line_in)
#,x,y2,x,y3,x,y4,x,y5,x,y6,x,y7,x,y8,x,y9,x,y10,'r')
	plt.title("Sensitivity of Gas Species to parameter "+k)
	ax3.legend(custom_lines,molecules,title="Gas Species")
	ax3.annotate("", xy=(0.5, 0.5), xytext=(0, 0),arrowprops=dict(arrowstyle="->"))
	ax3.set_xlabel('t (s)')
	ax3.set_ylabel(r'$\frac{ \partial{F} }{ \partial{k} } $',fontsize=15)

	parameters = np.linspace(1,10,10)
	norm = m_plt.colors.Normalize(vmin=np.min(parameters),vmax=np.max(parameters))
	s_m = m_plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	s_m.set_array([])

	divider = make_axes_locatable(ax3)
	cax = divider.append_axes("right", size="2%", pad=0.05)

	plt.colorbar(s_m,label='Pulse Number',cax=cax)
	plt.savefig('./'+output_file_name+'.png')
	#plt.show()