import numpy as np
import pandas as pd
import matplotlib as m_plt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

cmap = plt.get_cmap('rainbow')

custom_lines = [Line2D([0],[0],linestyle='--',color='b'), Line2D([0],[0],linestyle='-',color='orange'), Line2D([0],[0],linestyle='-.',color='g')]

constants = ['Ke0','Kd0','Ke1','Kd1','Ke2']
molecules = ['CO','O2','CO2']
#pulse_intensity = 2*1.6727777777777779e+19
folder_location = './'

x_location = [0.63,0.85,0.8,0.58,0.84]
y_location = [0.12,0.12,0.97,0.12,0.12]

ks = [1e-3,1e-4,1e-7,1,1]
#ks = [1,1,1,1,1]
elem_reactions = ['CO + * -> CO*','CO* -> CO + *','O2 + 2* -> 2O*','2O* -> O2 + 2*', 'CO* + O* -> CO2 + 2*']

df_co = pd.read_csv('../flux_data/CO.csv',header=None)
df_o2 = pd.read_csv('../flux_data/O2.csv',header=None)
df_co2 = pd.read_csv('../flux_data/CO2.csv',header=None)

fluxs = [df_co,df_o2,df_co2]

for k_num,k in enumerate(constants): #Generate a plot for each of the rate constants
	if k_num == 0:
		rcParams['figure.figsize'] = 2.75, 5
	else:
		rcParams['figure.figsize'] = 2, 5
	print(k)
	fig3, ax3 = plt.subplots()
	output_file_name = k+'_sens'
	for j_num,j in enumerate(molecules): # Add the lines for each one of the gaseous molecules
		print(j_num)
		gas = j
		parameter = k
		#Do the process for each of the pulses
		df_p1 = pd.read_csv(folder_location+gas+'/pulse_1.csv')

		
		x = df_p1['# t']
		
		#print(parameter)
		#print(df_p1['Ke0'])
		#print(df_p1[parameter])

		#sys.exit()
		y1 = df_p1[parameter]
		#print(y1)
		#y1*=ks[k_num]

		if j == 'CO':
			line_in = '--'
		elif j == 'O2':
			line_in = '-'
		else:
			line_in = '-.'
		#print((fluxs[j_num].iloc[1:,1]))
		#test = y1/(fluxs[j_num].iloc[0:,1])
		#print(y1)
		#sys.exit()
		#input_list = [y1/(fluxs[j_num].iloc[0:,1]),y2/(fluxs[j_num].iloc[1:,1]),y3/(fluxs[j_num].iloc[1:,1]),y4/(fluxs[j_num].iloc[1:,1]),y5/(fluxs[j_num].iloc[1:,1]),y6/(fluxs[j_num].iloc[1:,1]),y7/(fluxs[j_num].iloc[1:,1]),y8/(fluxs[j_num].iloc[1:,1]),y9/(fluxs[j_num].iloc[1:,1]),y10/(fluxs[j_num].iloc[1:,1])]

		#print(type(fluxs[j_num].iloc[1:,1]))
		#print(fluxs[j_num])
		#input_list = [y1/(fluxs[j_num].iloc[:,1])]
		input_list = [y1*(fluxs[j_num].iloc[:,1])/(fluxs[j_num].iloc[:,1])]
		
		#input_list = [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10]

		for k_du in range(0,len(input_list)):
			#print(input_list[k].iloc[0])
			input_list[k_du].iloc[0] = 0

		#print(input_list[0])
		#sys.exit()
		

		for z in range(0,1):
			ax3.plot(x,input_list[z],linestyle=line_in)
	#,x,y2,x,y3,x,y4,x,y5,x,y6,x,y7,x,y8,x,y9,x,y10,'r')
	#plt.title("Sensitivity of Gas Species to parameter "+k)
	props = dict(boxstyle='round', facecolor='white', alpha=0.4)
	#ax3.text(x_location[k_num], y_location[k_num], '\n'.join(('Constant: '+constants[k_num],elem_reactions[k_num])), transform=ax3.transAxes, fontsize=12,verticalalignment='top', bbox=props,ha='center')
	ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	#ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	#ax3.legend(custom_lines,molecules,title="Gas Species")
	ax3.annotate("", xy=(0.5, 0.5), xytext=(0, 0),arrowprops=dict(arrowstyle="->"))
	ax3.set_xlabel('t (s)')
	if k_num == 0:
		ax3.set_ylabel(r'$\frac{ \partial{F} }{ \partial{k} }$',fontsize=15)
		#ax3.set_ylabel(r'$ \frac{k}{F} (\frac{ \partial{F} }{ \partial{k} } )$',fontsize=15)

	parameters = np.linspace(1,10,10)
	norm = m_plt.colors.Normalize(vmin=np.min(parameters),vmax=np.max(parameters))
	#s_m = m_plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	#s_m.set_array([])
	ax3.set_xlim(0,0.07)
	divider = make_axes_locatable(ax3)
	#cax = divider.append_axes("right", size="2%", pad=0.05)
	plt.savefig('./'+k+'_sens.png')
	plt.clf()
	plt.close()
	#plt.show()
	#sys.exit()
	