import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import imageio
import pandas as pd
import sys

#rcParams.update({'figure.autolayout': True})
cmap = plt.get_cmap('Greys')
#custom_lines = [Line2D([0],[0],linestyle='--',color='k'), Line2D([0],[0],linestyle='-',color='k')]

exper_fol_loc = '../analysis/lin_fit_analysis/linear_model_folder/'

#exper_fol_loc = './paper_LH_exp_folder_folder/flux_data/'

#syn_fol_loc = '../analysis/CO_O2_analysis/CO_O2_fit_folder/fitting/iter_'

#syn_fol_loc = '../analysis/CO_O2_analysis/CO_O2_fit_3points_folder/fitting/iter_'

syn_fol_loc = '../analysis/lin_fit_analysis/gif_'

#syn_fol_loc = './LH_to_LH_1p_folder/fitting/iter_'

figure_name = './linear.png'

#figure_name = './ER_to_LH_1p.png'
#figure_name = './LH_to_LH.png' 

### FROM FILE

folders = [0,15,38]


#with open('../analysis/CO_O2_analysis/CO_O2_fit_folder/fitting/your_file.txt', 'r') as f:
#	lines = f.readlines()
#f.close
#lines = [x.strip() for x in lines]
#times = lines[0]
#times = times.replace('Contents: ','')
#
#times = eval(times)
#constants = lines[1]
#constants = constants.replace('Constants: ','')
#constants = eval(constants)
##things = len(times)

#first_it = constants[folders[0]]
#second_it = constants[folders[1]]
#third_it = constants[folders[2]]

f, axarr = plt.subplots(1, 3)

#print(constants[0])

#sys.exit()

### DIRECT

folders = [0,2,16]



first_it = [0.998,2.02]
#first_it = constants[folders[0]]#[0.0,0.0,0.0]
second_it = [0.998,2.02]
#second_it = constants[folders[1]]
#third_it =constants[folders[2]]

#first_it = [0,0,0]
#second_it = ['1.28e-3', '1.28e-3', '1.28e-3']
#third_it = [0.001,10,0.001]

third_it = [0.998,2.02]
#first_it = [1.71,1.71,1.71]
#second_it = [3.02,3.31,1.71]
#third_it = [8.68, 5.60, 1.71]

f.set_size_inches(8, 2)

#fig, ax = plt.subplots(figsize=(10,5))

props = dict(boxstyle='round', facecolor='white', alpha=0.4)

test = axarr[0]
test.set_title('Iteration: '+str(folders[0]))
test.set_ylabel('Flux')
test.set_xlabel('Time (s)')
folder = 0

#textstr = 'Ke0: '+str(second_it[0])
textstr = 'Ke0: '+str('{:.2e}'.format(first_it[0]))	
#textstr = '\n'.join((textstr,'D(Inert): '+str(first_it[0])))
#textstr = '\n'.join((textstr,'Kd0: '+str('{:.2e}'.format(first_it[1]))))
textstr = '\n'.join((textstr,'Ke1: '+str('{:.2e}'.format(first_it[1]))))
#textstr = '\n'.join((textstr,'Kd1: '+str('{:.2e}'.format(first_it[3]))))
#textstr = '\n'.join((textstr,'Ke2: '+str('{:.2e}'.format(first_it[4]))))


step_data_0 = pd.read_csv(syn_fol_loc+str(folders[0])+'_folder/flux_data/Inert.csv',header=None)
step_data_1 = pd.read_csv(syn_fol_loc+str(folders[0])+'_folder/flux_data/A.csv',header=None)
step_data_2 = pd.read_csv(syn_fol_loc+str(folders[0])+'_folder/flux_data/B.csv',header=None)
#step_data_3 = pd.read_csv(syn_fol_loc+str(folders[0])+'_folder/flux_data/O2.csv',header=None)

test.plot(step_data_0[0], step_data_0[1],ls='--',c='r')#,c=color
test.plot(step_data_1[0], step_data_1[1],ls='--',c='c')
test.plot(step_data_2[0], step_data_2[1],ls='--',c='b')
#test.plot(step_data_3[0], step_data_3[1],ls='--',c='g')

test.text(0.5, 0.87, textstr, transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')

test = axarr[1]
folder = 2
test.set_title('Iteration: '+str(folders[1]))
test.set_xlabel('Time (s)')

#textstr = 'Ke0: '+str(second_it[0])
textstr = 'Ke0: '+str('{:.2e}'.format(second_it[0]))	
#textstr = '\n'.join((textstr,'D(Inert): '+str(first_it[0])))
#textstr = '\n'.join((textstr,'Kd0: '+str('{:.2e}'.format(first_it[1]))))
textstr = '\n'.join((textstr,'Ke1: '+str('{:.2e}'.format(second_it[1]))))
#textstr = '\n'.join((textstr,'Kd1: '+str('{:.2e}'.format(first_it[3]))))
#textstr = '\n'.join((textstr,'Ke2: '+str('{:.2e}'.format(first_it[4]))))

step_data_0 = pd.read_csv(syn_fol_loc+str(folders[1])+'_folder/flux_data/Inert.csv',header=None)
step_data_1 = pd.read_csv(syn_fol_loc+str(folders[1])+'_folder/flux_data/A.csv',header=None)
step_data_2 = pd.read_csv(syn_fol_loc+str(folders[1])+'_folder/flux_data/B.csv',header=None)
#step_data_3 = pd.read_csv(syn_fol_loc+str(folders[1])+'_folder/flux_data/O2.csv',header=None)

test.plot(step_data_0[0], step_data_0[1],ls='--',c='r')#,c=color
test.plot(step_data_1[0], step_data_1[1],ls='--',c='c')
test.plot(step_data_2[0], step_data_2[1],ls='--',c='b')
#test.plot(step_data_3[0], step_data_3[1],ls='--',c='g')

test.text(0.5, 0.87, textstr, transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')


test = axarr[2]
folder = 6
test.set_title('Iteration: '+str(folders[2]))
test.set_xlabel('Time (s)')

#textstr = 'Ke0: '+str(second_it[0])
textstr = 'Ke0: '+str('{:.2e}'.format(third_it[0]))	
#textstr = '\n'.join((textstr,'D(Inert): '+str(first_it[0])))
#textstr = '\n'.join((textstr,'Kd0: '+str('{:.2e}'.format(first_it[1]))))
textstr = '\n'.join((textstr,'Ke1: '+str('{:.2e}'.format(third_it[1]))))
#textstr = '\n'.join((textstr,'Kd1: '+str('{:.2e}'.format(first_it[3]))))
#textstr = '\n'.join((textstr,'Ke2: '+str('{:.2e}'.format(first_it[4]))))


step_data_0 = pd.read_csv(syn_fol_loc+str(folders[1])+'_folder/flux_data/Inert.csv',header=None)
step_data_1 = pd.read_csv(syn_fol_loc+str(folders[1])+'_folder/flux_data/A.csv',header=None)
step_data_2 = pd.read_csv(syn_fol_loc+str(folders[1])+'_folder/flux_data/B.csv',header=None)
#step_data_3 = pd.read_csv(syn_fol_loc+str(folders[1])+'_folder/flux_data/O2.csv',header=None)

test.plot(step_data_0[0], step_data_0[1],ls='--',c='r')#,c=color
test.plot(step_data_1[0], step_data_1[1],ls='--',c='c')
test.plot(step_data_2[0], step_data_2[1],ls='--',c='b')
#test.plot(step_data_3[0], step_data_3[1],ls='--',c='g')

test.text(0.5, 0.87, textstr, transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')

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

exp_data_inert = pd.read_csv(exper_fol_loc+'Inert.csv',header=None)
exp_data_1 = pd.read_csv(exper_fol_loc+'A.csv',header=None)
exp_data_2 = pd.read_csv(exper_fol_loc+'B.csv',header=None)
#exp_data_3 = pd.read_csv(exper_fol_loc+'O2.csv',header=None)

test = axarr[0]
test.plot(exp_data_inert[0], exp_data_inert[1],c='r')
test.plot(exp_data_1[0], exp_data_1[1],c='c')
test.plot(exp_data_2[0], exp_data_2[1],c='b')
#test.plot(exp_data_3[0], exp_data_3[1],c='g')

#peak_loc = exp_data_inert.iloc[exp_data_inert[1].idxmax()]
#peak2 = exp_data_inert.loc[exp_data_inert[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')
#peak_loc = exp_data_1.iloc[exp_data_1[1].idxmax()]
#peak2 = exp_data_2.loc[exp_data_1[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')
#peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
#peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')



test = axarr[1]
test.plot(exp_data_inert[0], exp_data_inert[1],c='r')
test.plot(exp_data_1[0], exp_data_1[1],c='c')
test.plot(exp_data_2[0], exp_data_2[1],c='b')
#test.plot(exp_data_3[0], exp_data_3[1],c='g')

#peak_loc = exp_data_inert.iloc[exp_data_inert[1].idxmax()]
#peak2 = exp_data_inert.loc[exp_data_inert[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')
#peak_loc = exp_data_1.iloc[exp_data_1[1].idxmax()]
#peak2 = exp_data_2.loc[exp_data_1[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')
#peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
#peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')



test = axarr[2]
test.plot(exp_data_inert[0], exp_data_inert[1],c='r')
test.plot(exp_data_1[0], exp_data_1[1],c='c')
test.plot(exp_data_2[0], exp_data_2[1],c='b')
#test.plot(exp_data_3[0], exp_data_3[1],c='g')
#peak_loc = exp_data_inert.iloc[exp_data_inert[1].idxmax()]
#peak2 = exp_data_inert.loc[exp_data_inert[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')
#peak_loc = exp_data_1.iloc[exp_data_1[1].idxmax()]
#peak2 = exp_data_2.loc[exp_data_1[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')
#peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
#peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
#test.plot(peak_loc[0], peak_loc[1], 'ro')

f.tight_layout()	

plt.savefig(figure_name)
plt.show()