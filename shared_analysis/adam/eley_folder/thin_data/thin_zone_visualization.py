import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp
import imageio
import sys
import os

sim_info = pd.read_csv('../input_file.csv',header = None)
curr_time_steps = 250

length = float(sim_info[2][3])
mesh_size = float(sim_info[2][4])
time_steps = float(sim_info[2][2])
time_tot = float(sim_info[2][1])
cat_frac = float(sim_info[2][5])
reac_radius = float(sim_info[2][6])

species = 'O*'

thin_data = pd.read_csv('./'+species+'.csv',header = None)
#print(thin_data.iloc[2].sum()*)

x_length = mp.ceil(cat_frac*mesh_size)

x_dim = []

time = []
av_conc = []

for k in range(x_length-1):
	x_dim.append((k+1)*length/mesh_size)


def add_subplot_axes(ax,rect,axisbg='w'):
	fig = plt.gcf()
	box = ax.get_position()
	width = box.width
	height = box.height
	inax_position  = ax.transAxes.transform(rect[0:2])
	transFigure = fig.transFigure.inverted()
	infig_position = transFigure.transform(inax_position)    
	x = infig_position[0]
	y = infig_position[1]
	width *= rect[2]
	height *= rect[3]
	subax = fig.add_axes([x,y,width,height])
	x_labelsize = subax.get_xticklabels()[0].get_size()
	y_labelsize = subax.get_yticklabels()[0].get_size()
	x_labelsize *= rect[2]**0.5
	y_labelsize *= rect[3]**0.5
	subax.xaxis.set_tick_params(labelsize=x_labelsize)
	subax.yaxis.set_tick_params(labelsize=y_labelsize)
	return subax


def tap_plot(step):
	fig, ax = plt.subplots(figsize=(10,5))

	av_value = thin_data.iloc[step].sum()*(length/mesh_size)/(length*cat_frac)
	av_conc.append(av_value)
	time.append(step*(time_tot/time_steps))
	ax.grid()
	props = dict(facecolor='white')

	subpos = [0.4,0.35,0.5,0.5]
	#subpos = [0.2,0.6,0.3,0.3]
	subax1 = add_subplot_axes(ax,subpos)
	subax1.plot(time,av_conc)
	subax1.set_title('$Average\ TZ Value$')
	subax1.set_ylabel('$nmol/cm3$')
	subax1.set_xlabel('$t\ (seconds)$')
	subax1.set_xlim(0,0.15)
	subax1.set_ylim(0,10)

	ax.scatter(x_dim, thin_data.iloc[step])
	ax.set_xlabel('$Catalyst\ Zone\ (cm)$', fontsize=16)
	ax.set_ylabel('$Concentration\ (nmol/cm3)$', fontsize=16)
	
	ax.text(0.02, 0.95,'Time: '+str(round(step*(time_tot/time_steps),3))+' s',transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props)
	print(av_value)
	ax.hlines(y=av_value,xmin=(length/mesh_size),xmax=(x_length-1)*(length/mesh_size),linestyle='--',colors = 'b')
	#ax.set_ylim()
	ax.set_xlim(0, x_length*(length/mesh_size))
	ax.set_ylim(0, 20)
	#ax.set_ylim(270000, 280000)

	fig.canvas.draw()       # draw the canvas, cache the renderer
	image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
	image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

	return image

kwargs_write = {'fps':20.0, 'quantizer':'nq'}

imageio.mimsave('./'+species+'.gif', [tap_plot(i) for i in range(curr_time_steps)], fps=20)
