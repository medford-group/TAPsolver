import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import imageio
import pandas as pd
import sys
import random

things = [[0, 0, 0],
[1.31217914e-03, 2.63921671e-10, 6.10358601e-05],
[1.14680706e-03, 3.80401670e-10, 6.15471464e-05],
[1.20456432e-03, 3.03009877e-09, 6.86623943e-05],
[1.40565425e-03, 3.05425560e-08, 1.45134612e-04],
[1.32858346e-03, 4.40451496e-08, 1.83921683e-04],
[1.14331835e-03, 6.70134497e-08, 2.50403541e-04],
[1.11481281e-03, 8.90021543e-08, 3.12508878e-04],
[1.11717736e-03, 1.22408854e-07, 4.06153238e-04],
[1.08224268e-03, 1.62278844e-07, 5.18510638e-04],
[1.07411335e-03, 2.13139561e-07, 6.60066842e-04],
[1.02341593e-03, 2.32184121e-07, 7.34223943e-04],
[1.07216075e-03, 3.72658649e-07, 1.10900755e-03],
[1.04481282e-03, 4.23426579e-07, 1.26285025e-03],
[1.02286347e-03, 5.17281318e-07, 1.53549608e-03],
[1.01900882e-03, 6.14688858e-07, 1.81052488e-03],
[1.02088013e-03, 7.56012889e-07, 2.20264928e-03],
[1.02021746e-03, 9.14628973e-07, 2.64403097e-03],
[1.01698071e-03, 1.11711233e-06, 3.20871615e-03],
[1.01370826e-03, 1.36508208e-06, 3.89258407e-03],
[1.01150480e-03, 1.66638261e-06, 4.72016123e-03],
[1.01006777e-03, 2.05674810e-06, 5.72532919e-03],
[1.01003091e-03, 2.97597226e-06, 8.10522587e-03],
[1.01003091e-03, 2.97597226e-06, 8.10522587e-03],
[1.00638565e-03, 3.69008437e-06, 9.95464045e-03],
[1.00547267e-03, 4.37915554e-06, 1.17389536e-02],
[1.00642518e-03, 5.40837453e-06, 1.44074455e-02],
[1.00589342e-03, 6.51744850e-06, 1.72823979e-02],
[1.00553988e-03, 6.53161721e-06, 1.73176871e-02],
[1.00505870e-03, 8.09618791e-06, 2.13738904e-02],
[1.00465838e-03, 9.92164469e-06, 2.61066139e-02],
[1.00437594e-03, 1.22249053e-05, 3.20787859e-02],
[1.00414553e-03, 1.50588742e-05, 3.94274883e-02],
[1.00395600e-03, 1.86044057e-05, 4.86216936e-02],
[1.00380107e-03, 2.30362356e-05, 6.01145277e-02],
[1.00367428e-03, 2.86120643e-05, 7.45742694e-02],
[1.00357084e-03, 3.56566371e-05, 9.28430831e-02],
[1.00348664e-03, 4.46062767e-05, 1.16052488e-01],
[1.00341835e-03, 5.60370381e-05, 1.45696394e-01],
[1.00345242e-03, 5.60210807e-05, 1.45655012e-01],
[1.00337824e-03, 7.07521314e-05, 1.83857705e-01],
[1.00331553e-03, 8.98704489e-05, 2.33438095e-01],
[1.00328088e-03, 1.14507199e-04, 2.97329684e-01],
[1.00325342e-03, 1.46623859e-04, 3.80619261e-01],
[1.00323167e-03, 1.88659511e-04, 4.89632201e-01],
[1.00321467e-03, 2.43856382e-04, 6.32776724e-01],
[1.00320148e-03, 3.16523227e-04, 8.21226913e-01],
[1.00319130e-03, 4.12378392e-04, 1.06981241e+00],
[1.00318348e-03, 5.39001896e-04, 1.39819082e+00],
[1.00317751e-03, 7.06436926e-04, 1.83240759e+00]]


#sys.exit()
#def generate_gif(exp_data, iterations,molecules,peaks,equations,constants):
def generate_gif(all_steps):
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
		height *= rect[3]  # <= Typo was here
		subax = fig.add_axes([x,y,width,height])
		x_labelsize = subax.get_xticklabels()[0].get_size()
		y_labelsize = subax.get_yticklabels()[0].get_size()
		x_labelsize *= rect[2]**0.5
		y_labelsize *= rect[3]**0.5
		subax.xaxis.set_tick_params(labelsize=x_labelsize)
		subax.yaxis.set_tick_params(labelsize=y_labelsize)
		return subax
	x_data = list(range(0, all_steps))
	y_data = x_data

	def tap_plot(step):
		fig, ax = plt.subplots(figsize=(10,5))
	
		ax.grid()
		ax.set_ylim(0,0.014)
		#exp_data_inert = pd.read_csv('./paper_LH_exp_folder/flux_data/Inert.csv',header=None)
		exp_data_1 = pd.read_csv('./paper_ER_exp_folder/flux_data/CO'+'.csv',header=None)
		exp_data_2 = pd.read_csv('./paper_ER_exp_folder/flux_data/CO2'+'.csv',header=None)
		
		#ax.plot(exp_data_inert[0], exp_data_inert[1])
		ax.plot(exp_data_1[0], exp_data_1[1])
		ax.plot(exp_data_2[0], exp_data_2[1])
	
		#step_data_1 = pd.read_csv('./gif_'+str(step)+'_folder/flux_data/A'+'.csv',header=None)
		#step_data_2 = pd.read_csv('./gif_'+str(step)+'_folder/flux_data/B'+'.csv',header=None)
	
		#step_data_0 = pd.read_csv('./diff_data_folder/flux_data/Inert'+'.csv',header=None)
		#step_data_0 = pd.read_csv('./multi_mechanisms/lh_to_lh_'+str(step)+'_folder/flux_data/Inert'+'.csv',header=None)
		#multi_mechanisms/

		step_data_1 = pd.read_csv('./five_point_lh_to_eh_'+str(step)+'_folder/flux_data/CO'+'.csv',header=None)
		step_data_2 = pd.read_csv('./five_point_lh_to_eh_'+str(step)+'_folder/flux_data/CO2'+'.csv',header=None)
	
		#step_data_1 = pd.read_csv('./diff_data_folder/flux_data/A'+'.csv',header=None)
		#step_data_2 = pd.read_csv('./diff_data_folder/flux_data/B'+'.csv',header=None)
	
		#ax.plot(step_data_0[0], step_data_0[1],ls='--')
		ax.plot(step_data_1[0], step_data_1[1],ls='--')
		ax.plot(step_data_2[0], step_data_2[1],ls='--')
		ax.set_xlabel('Time (s)', fontsize=16)
		ax.set_ylabel('Flux (nmol/s)', fontsize=16)
		textstr = 'Rate Constants: '	
		#textstr = '\n'.join((textstr,'D(Inert): '+str(first_it[0])))
		textstr = '\n'.join((textstr,'kf1: '+'{:0.3e}'.format(things[step][0])))
		textstr = '\n'.join((textstr,'kb1: '+'{:0.3e}'.format(things[step][1])))
		textstr = '\n'.join((textstr,'kf2: '+'{:0.3e}'.format(things[step][2])))
		props = dict(facecolor='white')
		ax.text(0.8, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props,ha='center')
		textstr = 'Elementary Reactions:'	
		#textstr = '\n'.join((textstr,'D(Inert): '+str(first_it[0])))
		textstr = '\n'.join((textstr,'CO + * <-> CO*'))
		textstr = '\n'.join((textstr,'CO* + O* -> CO2 + 2*'))
		ax.text(0.5, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props,ha='center')
		ax.text(0.15, 0.95,'Iteration: '+str(step),transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props)
		#peak_loc = exp_data_inert.iloc[exp_data_inert[1].idxmax()]
		#peak2 = exp_data_inert.loc[exp_data_inert[0] == peak_loc[0]].index
		#plt.plot(peak_loc[0], peak_loc[1], 'ro')
		peak_loc = exp_data_1.iloc[exp_data_1[1].idxmax()]
		peak2 = exp_data_2.loc[exp_data_1[0] == peak_loc[0]].index
		plt.plot(peak_loc[0], peak_loc[1], 'ro')
		peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
		peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
		plt.plot(peak_loc[0], peak_loc[1], 'ro')
		#print(peak_loc)
		#print(peak2)
		#sys.exit()
		#x = np.linspace(-np.pi,np.pi)
		ax.set_ylim(0,0.014)
		subpos = [0.4,0.13,0.5,0.4]
		#subpos = [0.2,0.6,0.3,0.3]
		subax1 = add_subplot_axes(ax,subpos)
		subax1.plot(x_data[:step],y_data[:step])
		subax1.set_title('Time per Iteration')
		subax1.set_ylabel('Time (minutes)')
		subax1.set_xlabel('Iteration #')
		subax1.set_xlim(0,all_steps)
		subax1.set_ylim(0,all_steps)
		#ax.set_ylim(0, y_max)
		fig.canvas.draw()       # draw the canvas, cache the renderer
		image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
		image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	
		return image
	
		#plt.show()
	
	
	kwargs_write = {'fps':4.0, 'quantizer':'nq'}
	imageio.mimsave('./test.gif', [tap_plot(i) for i in range(all_steps)], fps=4)
tot_steps = 49
#generate_gif(tot_steps)

#tap_plot(1)

def plot_for_offset(power, y_max):
	# Data for plotting
	t = np.arange(0.0, 100, 1)
	s = t**power

	fig, ax = plt.subplots(figsize=(10,5))
	ax.plot(t, s)
	ax.grid()
	ax.set(xlabel='X', ylabel='x',
		   title='Powers of x')

	# IMPORTANT ANIMATION CODE HERE
	# Used to keep the limits constant
	ax.set_ylim(0, y_max)

	# Used to return the plot as an image rray
	fig.canvas.draw()       # draw the canvas, cache the renderer
	image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
	image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

	return image

#kwargs_write = {'fps':20.0, 'quantizer':'nq'}
#imageio.mimsave('./powers.gif', [plot_for_offset(i, 100) for i in range(100)], fps=20)