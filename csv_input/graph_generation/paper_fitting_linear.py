import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import imageio
import pandas as pd

#rcParams.update({'figure.autolayout': True})
cmap = plt.get_cmap('Greys')
#custom_lines = [Line2D([0],[0],linestyle='--',color='k'), Line2D([0],[0],linestyle='-',color='k')]

f, axarr = plt.subplots(1, 3)

first_it = [0.5,0.77,1.71]
second_it = [0.89,1.4,1.71]
third_it = [0.99, 2.01, 1.71]

f.set_size_inches(8, 2)

#fig, ax = plt.subplots(figsize=(10,5))

props = dict(boxstyle='round', facecolor='white', alpha=0.4)

test = axarr[0]
test.set_title('Initial Guess')
test.set_ylabel('Flux')
test.set_xlabel('Time (s)')
folder = 0

textstr = 'k1: '+str(first_it[0])	
#textstr = '\n'.join((textstr,'D(Inert): '+str(first_it[0])))
textstr = '\n'.join((textstr,'k2: '+str(first_it[1])))

step_data_0 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/Inert.csv',header=None)
step_data_1 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/A.csv',header=None)
step_data_2 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/B.csv',header=None)

test.plot(step_data_0[0], step_data_0[1],ls='--',c='r')#,c=color
test.plot(step_data_1[0], step_data_1[1],ls='--',c='c')
test.plot(step_data_2[0], step_data_2[1],ls='--',c='b')

test.text(0.8, 0.87, textstr, transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')

test = axarr[1]
folder = 10
test.set_title('Iteration: '+str(folder))
test.set_xlabel('Time (s)')

textstr = 'k1: '+str(second_it[0])	
#textstr = '\n'.join((textstr,'D(Inert): '+str(first_it[0])))
textstr = '\n'.join((textstr,'k2: '+str(second_it[1])))

step_data_0 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/Inert.csv',header=None)
step_data_1 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/A.csv',header=None)
step_data_2 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/B.csv',header=None)

test.plot(step_data_0[0], step_data_0[1],ls='--',c='r')#,c=color
test.plot(step_data_1[0], step_data_1[1],ls='--',c='c')
test.plot(step_data_2[0], step_data_2[1],ls='--',c='b')

test.text(0.8, 0.87, textstr, transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')


test = axarr[2]
folder = 15
test.set_title('Iteration: '+str(folder))
test.set_xlabel('Time (s)')

textstr = 'k1: '+str(third_it[0])	
#textstr = '\n'.join((textstr,'D(Inert): '+str(first_it[0])))
textstr = '\n'.join((textstr,'k2: '+str(third_it[1])))

step_data_0 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/Inert.csv',header=None)
step_data_1 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/A.csv',header=None)
step_data_2 = pd.read_csv('../analysis/lin_fit_analysis/gif_'+str(folder)+'_folder/flux_data/B.csv',header=None)

test.plot(step_data_0[0], step_data_0[1],ls='--',c='r')#,c=color
test.plot(step_data_1[0], step_data_1[1],ls='--',c='c')
test.plot(step_data_2[0], step_data_2[1],ls='--',c='b')

test.text(0.8, 0.87, textstr, transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')

#for step_2 in range(0, 1):
	
#	step = step_2
	#step_data_1 = pd.read_csv('./gif_'+str(step)+'_folder/flux_data/A'+'.csv',header=None)
	#step_data_2 = pd.read_csv('./gif_'+str(step)+'_folder/flux_data/B'+'.csv',header=None)

	#step_data_0 = pd.read_csv('./diff_data_folder/flux_data/Inert'+'.csv',header=None)
#	step_data_0 = pd.read_csv('./alt_units_fit_folder/flux_data/Inert.csv',header=None)

#	step_data_1 = pd.read_csv('./alt_units_fit_folder/flux_data/CO.csv',header=None)
#	step_data_2 = pd.read_csv('./alt_units_fit_folder/flux_data/CO2.csv',header=None)

	#step_data_1 = pd.read_csv('./diff_data_folder/flux_data/A'+'.csv',header=None)
	#step_data_2 = pd.read_csv('./diff_data_folder/flux_data/B'+'.csv',header=None)
	#color = cmap(float(step)/15)

#	ax.plot(step_data_0[0], step_data_0[1],ls='--')#,c=color
#	ax.plot(step_data_1[0], step_data_1[1],ls='--')
#	ax.plot(step_data_2[0], step_data_2[1],ls='--')


	#print(peak_loc)
	#print(peak2)
	#sys.exit()
	
	#ax.set_ylim(0, y_max)
	#fig.canvas.draw()       # draw the canvas, cache the renderer
	#image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
	#image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

	

	#plt.show()

exp_data_inert = pd.read_csv('../analysis/lin_fit_analysis/gif_18_folder/flux_data/Inert.csv',header=None)
exp_data_1 = pd.read_csv('../analysis/lin_fit_analysis/gif_18_folder/flux_data/A'+'.csv',header=None)
exp_data_2 = pd.read_csv('../analysis/lin_fit_analysis/gif_18_folder/flux_data/B'+'.csv',header=None)

test = axarr[0]
test.plot(exp_data_inert[0], exp_data_inert[1],c='r')
test.plot(exp_data_1[0], exp_data_1[1],c='c')
test.plot(exp_data_2[0], exp_data_2[1],c='b')


peak_loc = exp_data_1.iloc[exp_data_1[1].idxmax()]
peak2 = exp_data_2.loc[exp_data_1[0] == peak_loc[0]].index
test.plot(peak_loc[0], peak_loc[1], 'ro')
peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
test.plot(peak_loc[0], peak_loc[1], 'ro')



test = axarr[1]
test.plot(exp_data_inert[0], exp_data_inert[1],c='r')
test.plot(exp_data_1[0], exp_data_1[1],c='c')
test.plot(exp_data_2[0], exp_data_2[1],c='b')


peak_loc = exp_data_1.iloc[exp_data_1[1].idxmax()]
peak2 = exp_data_2.loc[exp_data_1[0] == peak_loc[0]].index
test.plot(peak_loc[0], peak_loc[1], 'ro')
peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
test.plot(peak_loc[0], peak_loc[1], 'ro')



test = axarr[2]
test.plot(exp_data_inert[0], exp_data_inert[1],c='r')
test.plot(exp_data_1[0], exp_data_1[1],c='c')
test.plot(exp_data_2[0], exp_data_2[1],c='b')

peak_loc = exp_data_1.iloc[exp_data_1[1].idxmax()]
peak2 = exp_data_2.loc[exp_data_1[0] == peak_loc[0]].index
test.plot(peak_loc[0], peak_loc[1], 'ro')
peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
test.plot(peak_loc[0], peak_loc[1], 'ro')

f.tight_layout()	
plt.savefig('./lin_paper.png')
plt.show()