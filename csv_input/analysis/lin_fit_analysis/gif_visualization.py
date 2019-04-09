import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import imageio
import pandas as pd


def tap_plot(step):
	fig, ax = plt.subplots(figsize=(10,5))
	
	ax.grid()

	exp_data_inert = pd.read_csv('./exp_data_folder/flux_data/Inert.csv',header=None)
	exp_data_1 = pd.read_csv('./exp_data_folder/flux_data/CO'+'.csv',header=None)
	exp_data_2 = pd.read_csv('./exp_data_folder/flux_data/CO2'+'.csv',header=None)
	
	ax.plot(exp_data_inert[0], exp_data_inert[1])
	ax.plot(exp_data_1[0], exp_data_1[1])
	ax.plot(exp_data_2[0], exp_data_2[1])

	#step_data_1 = pd.read_csv('./gif_'+str(step)+'_folder/flux_data/A'+'.csv',header=None)
	#step_data_2 = pd.read_csv('./gif_'+str(step)+'_folder/flux_data/B'+'.csv',header=None)

	#step_data_0 = pd.read_csv('./diff_data_folder/flux_data/Inert'+'.csv',header=None)
	step_data_0 = pd.read_csv('./lh_to_lh_'+str(step)+'_folder/flux_data/Inert'+'.csv',header=None)

	step_data_1 = pd.read_csv('./lh_to_lh_'+str(step)+'_folder/flux_data/CO'+'.csv',header=None)
	step_data_2 = pd.read_csv('./lh_to_lh_'+str(step)+'_folder/flux_data/CO2'+'.csv',header=None)

	#step_data_1 = pd.read_csv('./diff_data_folder/flux_data/A'+'.csv',header=None)
	#step_data_2 = pd.read_csv('./diff_data_folder/flux_data/B'+'.csv',header=None)

	ax.plot(step_data_0[0], step_data_0[1],ls='--')
	ax.plot(step_data_1[0], step_data_1[1],ls='--')
	ax.plot(step_data_2[0], step_data_2[1],ls='--')

	peak_loc = exp_data_inert.iloc[exp_data_inert[1].idxmax()]
	peak2 = exp_data_inert.loc[exp_data_inert[0] == peak_loc[0]].index
	plt.plot(peak_loc[0], peak_loc[1], 'ro')
	peak_loc = exp_data_1.iloc[exp_data_1[1].idxmax()]
	peak2 = exp_data_2.loc[exp_data_1[0] == peak_loc[0]].index
	plt.plot(peak_loc[0], peak_loc[1], 'ro')
	peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
	peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
	plt.plot(peak_loc[0], peak_loc[1], 'ro')
	#print(peak_loc)
	#print(peak2)
	#sys.exit()
	
	#ax.set_ylim(0, y_max)
	fig.canvas.draw()       # draw the canvas, cache the renderer
	image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
	image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

	return image

	#plt.show()


kwargs_write = {'fps':3.0, 'quantizer':'nq'}
imageio.mimsave('./lh_to_lh_tap.gif', [tap_plot(i) for i in range(26)], fps=3)

tap_plot(1)

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