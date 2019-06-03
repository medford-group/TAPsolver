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

custom_lines = [Line2D([0],[0],linestyle='--',color='k'), Line2D([0],[0],linestyle='-',color='k')]

constants = ['Ke0','Kd0','Ke1','Ke2','Ke3']
molecules = ['CO','CO2']
pulse_intensity = 2*1.6727777777777779e+19
folder_location = './'

x_location = [0.63,0.85,0.8,0.58,0.84]
y_location = [0.12,0.12,0.97,0.12,0.12]

ks = [1.1108e07,1.6900e12,3.354e-7,3.287e-10,9.6077e9]
#ks = [1,1,1,1,1]
elem_reactions = ['CO -> CO*','CO* -> CO','CO* + OA* -> CO2*','CO* + OB* -> CO2*', 'CO2* -> CO2']

df_co = pd.read_csv('../CO.csv',header=None)
df_co2 = pd.read_csv('../CO2.csv',header=None)

fluxs = [df_co,df_co2]

for k_num,k in enumerate(constants): #Generate a plot for each of the rate constants
	fig3, ax3 = plt.subplots()
	output_file_name = k+'_sens'
	for j_num,j in enumerate(molecules): # Add the lines for each one of the gaseous molecules
		gas = j
		parameter = k
		#Do the process for each of the pulses
		df_p1 = pd.read_csv(folder_location+gas+'/pulse_1.csv')
		df_p2 = pd.read_csv(folder_location+gas+'/pulse_2.csv')
		df_p3 = pd.read_csv(folder_location+gas+'/pulse_3.csv')
		df_p4 = pd.read_csv(folder_location+gas+'/pulse_4.csv')
		df_p5 = pd.read_csv(folder_location+gas+'/pulse_5.csv')
		df_p6 = pd.read_csv(folder_location+gas+'/pulse_6.csv')
		df_p7 = pd.read_csv(folder_location+gas+'/pulse_7.csv')
		df_p8 = pd.read_csv(folder_location+gas+'/pulse_8.csv')
		df_p9 = pd.read_csv(folder_location+gas+'/pulse_9.csv')
		df_p10 = pd.read_csv(folder_location+gas+'/pulse_10.csv')
		
		x = df_p1['# t']
		
		#print(parameter)
		#print(df_p1['Ke0'])
		#print(df_p1[parameter])

		#sys.exit()
		y1 = df_p1[parameter]/pulse_intensity
		print(y1)
		y1*=ks[k_num]
		
		#sys.exit()
		y2 = df_p2[parameter]/pulse_intensity
		y2*=ks[k_num]
		y3 = df_p3[parameter]/pulse_intensity
		y3*=ks[k_num]
		y4 = df_p4[parameter]/pulse_intensity
		y4*=ks[k_num]
		y5 = df_p5[parameter]/pulse_intensity
		y5*=ks[k_num]
		y6 = df_p6[parameter]/pulse_intensity
		y6*=ks[k_num]
		y7 = df_p7[parameter]/pulse_intensity
		y7*=ks[k_num]
		y8 = df_p8[parameter]/pulse_intensity
		y8*=ks[k_num]
		y9 = df_p9[parameter]/pulse_intensity
		y9*=ks[k_num]
		y10 = df_p10[parameter]/pulse_intensity
		y10*=ks[k_num]
		if j == 'CO':
			line_in = '--'
		else:
			line_in = '-'
		#print((fluxs[j_num].iloc[1:,1]))
		#test = y1/(fluxs[j_num].iloc[0:,1])
		#print(y1)
		#sys.exit()
		#input_list = [y1/(fluxs[j_num].iloc[0:,1]),y2/(fluxs[j_num].iloc[1:,1]),y3/(fluxs[j_num].iloc[1:,1]),y4/(fluxs[j_num].iloc[1:,1]),y5/(fluxs[j_num].iloc[1:,1]),y6/(fluxs[j_num].iloc[1:,1]),y7/(fluxs[j_num].iloc[1:,1]),y8/(fluxs[j_num].iloc[1:,1]),y9/(fluxs[j_num].iloc[1:,1]),y10/(fluxs[j_num].iloc[1:,1])]
		print(fluxs[0].iloc[:,1])
		print(fluxs[1].iloc[:,1])
		sys.exit()
		#print(type(fluxs[j_num].iloc[1:,1]))
		#print(fluxs[j_num])
		input_list = [y1/(fluxs[j_num].iloc[:,1]),y2/(fluxs[j_num].iloc[1:,2]),y3/(fluxs[j_num].iloc[1:,3]),y4/(fluxs[j_num].iloc[1:,4]),y5/(fluxs[j_num].iloc[1:,5]),y6/(fluxs[j_num].iloc[1:,6]),y7/(fluxs[j_num].iloc[1:,7]),y8/(fluxs[j_num].iloc[1:,8]),y9/(fluxs[j_num].iloc[1:,9]),y10/(fluxs[j_num].iloc[1:,10])]
		
		#input_list = [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10]

		for k_du in range(0,len(input_list)):
			#print(input_list[k].iloc[0])
			input_list[k_du].iloc[0] = 0

		#print(input_list[0])
		#sys.exit()
		

		for z in range(0,10):
			color = cmap(float(z)/10)
			ax3.plot(x,input_list[z],c=color,linestyle=line_in)
	#,x,y2,x,y3,x,y4,x,y5,x,y6,x,y7,x,y8,x,y9,x,y10,'r')
	#plt.title("Sensitivity of Gas Species to parameter "+k)
	props = dict(boxstyle='round', facecolor='white', alpha=0.4)
	ax3.text(x_location[k_num], y_location[k_num], '\n'.join(('Constant: '+constants[k_num],elem_reactions[k_num])), transform=ax3.transAxes, fontsize=12,verticalalignment='top', bbox=props,ha='center')
	ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax3.legend(custom_lines,molecules,title="Gas Species")
	ax3.annotate("", xy=(0.5, 0.5), xytext=(0, 0),arrowprops=dict(arrowstyle="->"))
	ax3.set_xlabel('t (s)')
	ax3.set_ylabel(r'$ \frac{k}{F} (\frac{ \partial{F} }{ \partial{k} } )$',fontsize=15)

	parameters = np.linspace(1,10,10)
	norm = m_plt.colors.Normalize(vmin=np.min(parameters),vmax=np.max(parameters))
	s_m = m_plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	s_m.set_array([])

	divider = make_axes_locatable(ax3)
	cax = divider.append_axes("right", size="2%", pad=0.05)

	plt.colorbar(s_m,label='Pulse Number',cax=cax)
	#plt.savefig('./'+output_file_name+'.png')
	plt.show()
	#sys.exit()