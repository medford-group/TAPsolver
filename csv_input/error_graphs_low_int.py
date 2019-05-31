import sys
import pandas as pd
#from user_input import base_run
import imageio
import matplotlib.pyplot as plt
import numpy as np

syn_data_folder = './analysis/CO_O2_analysis/error_v_ratio_low_int/'#'./analysis/CO_O2_analysis/CO_O2_exp_folder/flux_data/'
#two_syn_data_folder = './analysis/lh_eh_co_analysis/paper_ER_exp_folder/flux_data/'
exp_data_folder = './analysis/CO_O2_analysis/error_v_ratio_low_int/'#'./analysis/CO_O2_analysis/CO_O2_fit_productobj_3rd_folder/fitting/'
num_in_folder = 8
reactants = ['CO','O2','CO2','Inert']
x_data = [0.2/1.8,0.4/1.6,0.6/1.4,0.8/1.2,1,1.2/0.8,1.4/0.6,1.6/0.4,1.8/.2,]
y_data = [1.2220057569198642e-07,9.763829356572386e-08,7.741030285747366e-08,5.899042142877121e-08,4.331773234205131e-08,4.1852377644395655e-08,3.474962728673667e-08,2.537996426317133e-08,1.2584598587668078e-08]
def generate_gif(molecules,exp_loc,fit_loc,all_steps,constants,reactions,time_data,actual_mech):
	#constant_names = ['$k_{f1}$','$k_{b1}$','$k_{f2}$','$k_{b2}$','$k_{f3}$']
	#constant_names = ['M(A)','M(B)']
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
	for k in range(0,all_steps-1):
		print(k)
		#y_data.append( (time_data[k+1] - time_data[k]) / 60 )

	def tap_plot(step):
		
		colors = ['b','r','g','m','k','y','c']
		point_colors = ['bo','ro','go','mo','ko','yo','co']
		fig, ax = plt.subplots(figsize=(10,5))
	
		ax.grid(alpha=0.4)
		
		#exp_data_inert = pd.read_csv('./paper_LH_exp_folder/flux_data/Inert.csv',header=None)
		
		exp_data = {}
		#exp_data_2 = {}
		sim_data = {}
		print(step)
		for k_names in molecules:
			exp_data[k_names] = pd.read_csv(fit_loc+str(round(.2*(step+1),1))+'_'+str(round(2 - .2*(step+1),1))+'_folder/flux_data/'+k_names+'.csv',header=None)
			#exp_data_2[k_names] = pd.read_csv(two_syn_data_folder+k_names+'.csv',header=None)
			sim_data[k_names] = pd.read_csv(fit_loc+str(round(.2*(step+1),1))+'_'+str(round(2 - .2*(step+1),1))+'_folder/flux_data/'+k_names+'.csv',header=None)

		#exp_data_2 = pd.read_csv('./paper_ER_exp_folder/flux_data/CO2'+'.csv',header=None)
		
		#ax.plot(exp_data_inert[0], exp_data_inert[1])
		for k,k_names in enumerate(molecules):
			#newList1 = [x / 1e18 for x in exp_data[k_names][1]]
			#ax.plot(exp_data[k_names][0], newList1,c=colors[k], label=molecules[k])
			ax.plot(exp_data[k_names][0], exp_data[k_names][1],c=colors[k],label=molecules[k])
			#ax.plot(exp_data_2[k_names][0], exp_data_2[k_names][1],c=colors[k],label=molecules[k])
		#ax.plot(exp_data_2[0], exp_data_2[1])
		
		for k,k_names in enumerate(molecules):
			#newList1 = [x * NUMBER for x in sim_data[k_names][1]]
			#ax.plot(sim_data[k_names][0], newList1,ls='--',c=colors[k], label='_nolegend_')
			ax.plot(sim_data[k_names][0], sim_data[k_names][1],ls='--',c=colors[k], label='_nolegend_')
		#step_data_1 = pd.read_csv('./gif_'+str(step)+'_folder/flux_data/A'+'.csv',header=None)
		#step_data_2 = pd.read_csv('./gif_'+str(step)+'_folder/flux_data/B'+'.csv',header=None)
	
		#step_data_0 = pd.read_csv('./diff_data_folder/flux_data/Inert'+'.csv',header=None)
		#step_data_0 = pd.read_csv('./multi_mechanisms/lh_to_lh_'+str(step)+'_folder/flux_data/Inert'+'.csv',header=None)
		#multi_mechanisms/

		#step_data_1 = pd.read_csv('./five_point_lh_to_eh_'+str(step)+'_folder/flux_data/CO'+'.csv',header=None)
		#step_data_2 = pd.read_csv('./five_point_lh_to_eh_'+str(step)+'_folder/flux_data/CO2'+'.csv',header=None)
	
		#step_data_1 = pd.read_csv('./diff_data_folder/flux_data/A'+'.csv',header=None)
		#step_data_2 = pd.read_csv('./diff_data_folder/flux_data/B'+'.csv',header=None)
	
		#ax.plot(step_data_0[0], step_data_0[1],ls='--')
		#ax.plot(step_data_1[0], step_data_1[1],ls='--')
		#ax.plot(step_data_2[0], step_data_2[1],ls='--')
		ax.set_xlabel('Time (s)', fontsize=16)
		ax.set_ylabel('Flux (nmol/s)', fontsize=16)



		#peak_loc = exp_data_inert.iloc[exp_data_inert[1].idxmax()]
		#peak2 = exp_data_inert.loc[exp_data_inert[0] == peak_loc[0]].index
		#plt.plot(peak_loc[0], peak_loc[1], 'ro')
		peak_loc = exp_data[molecules[-1:][0]].iloc[exp_data[molecules[-1:][0]][1].idxmax()]
		peak_peak = peak_loc[1]
		#ax.legend(title="Gas Species", fontsize=14)
		ax.legend(fontsize=14,loc=2, bbox_to_anchor=(0.14, 0.35, 0.5, 0.5),framealpha=1)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		#for k,k_names in enumerate(molecules):
		for k,k_names in enumerate(molecules[:-1]):
			
			peak_loc = exp_data[k_names].iloc[exp_data[k_names][1].idxmax()]
			#peak2 = exp_data_2.loc[exp_data[k_names][0] == peak_loc[0]].index
			plt.plot(peak_loc[0], peak_loc[1], point_colors[k])
			
			#peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
			#peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
			#plt.plot(peak_loc[0], peak_loc[1], 'ro')
			#times.append(round(user_data[k_new].iloc[1,0],6))
			#values.append(user_data[k_new].iloc[1,1])

			peak2 = exp_data[k_names].loc[exp_data[k_names][0] == peak_loc[0]].index
			test3 = int(round((peak2[0]+1)/2,0))
			mid_loc = exp_data[k_names].iloc[test3,:]
			plt.plot(mid_loc[0], mid_loc[1], point_colors[k])

			plt.plot(round(exp_data[k_names].iloc[1,0],6), exp_data[k_names].iloc[1,1], point_colors[k])
			
			value_test = 0.85*peak_loc[1]
			sort_fif = exp_data[k_names].iloc[(exp_data[k_names][1]-value_test).abs().argsort()[:]]

			for k_j in range(0,sort_fif.index[0]):
				if sort_fif.iloc[k_j,0] > peak_loc[0]:
					fif_point = sort_fif.iloc[k_j]
					break
				else:
					pass

			test2 = exp_data[k_names].loc[exp_data[k_names][0] == fif_point[0]].index

			plt.plot(round(fif_point[0],6),fif_point[1], point_colors[k])
			#time_step.append(test2[0])
			#times.append(round(fif_point[0],6))
			#values.append(fif_point[1])
			#%20 point
			value_test = 0.2*peak_loc[1]
			sort_twen = exp_data[k_names].iloc[(exp_data[k_names][1]-value_test).abs().argsort()[:]]
			
			for k_j in range(0,sort_twen.index[0]):
				if sort_twen.iloc[k_j,0] > peak_loc[0]:
					twen_point = sort_twen.iloc[k_j]
					break
				else:
					pass

			test2 = exp_data[k_names].loc[exp_data[k_names][0] == twen_point[0]].index
			plt.plot(round(twen_point[0],6),twen_point[1], point_colors[k])
			#time_step.append(test2[0])
			#times.append(round(twen_point[0],6))
			#values.append(twen_point[1])

		subpos = [0.4,0.23,0.5,0.5]
		#subpos = [0.2,0.6,0.3,0.3]
		subax1 = add_subplot_axes(ax,subpos)
		#subax1.plot(x_data[:step],y_data[:step+1])
		subax1.scatter(x_data[:step],y_data[:step])
		#subax1.set_title('Time per Iteration')
		subax1.set_ylabel('Objective Value')
		subax1.set_xlabel('CO / O2 Feed Ratio')
		subax1.set_xlim(0,5)
		subax1.set_ylim(0,1.5e-07)


		ax.set_ylim(0,3e-2)
		#print(peak_loc)
		#print(peak2)
		#sys.exit()
		#x = np.linspace(-np.pi,np.pi)
		#ax.set_ylim(0,0.014)
		subpos = [0.4,0.13,0.5,0.4]
		#subpos = [0.2,0.6,0.3,0.3]
		##subax1 = add_subplot_axes(ax,subpos)
		##subax1.plot(x_data[:step],y_data[:step])
		##subax1.set_title('Time per Iteration')
		##subax1.set_ylabel('Time (minutes)')
		##subax1.set_xlabel('Iteration #')
		##subax1.set_xlim(0,all_steps)
		##subax1.set_ylim(0,max(y_data)*1.1)
		#ax.set_ylim(0, y_max)
		fig.canvas.draw()       # draw the canvas, cache the renderer
		image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
		image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	
		return image
	
		#plt.show()
	
	
	kwargs_write = {'fps':4.0, 'quantizer':'nq'}

	my_list = list(range(0, all_steps))
	for k in range(0,4):
		my_list.append(all_steps)
	imageio.mimsave(fit_loc+'/output_new.gif', [tap_plot(i) for i in my_list], fps=4)
	#imageio.mimsave(fit_loc+'/output_new.gif', [tap_plot(i) for i in range(all_steps)], fps=4)




#constants = [[600. ,  2.],[600.00714538,   2.99997447],[600.02740511  , 3.94722187],[600.05141097  , 4.2588502 ],[600.07777133 ,  4.32802698],[600.10630244 ,  4.33914975],[600.36227098 ,  4.38806658],[600.963561  ,   4.44768818],[602.85130388  , 4.55463656],[607.72506832 ,  4.71507406],[620.34409699 ,  4.95594844],[649.81348071,   5.26045016],[710.07998561 ,  5.51469226],[812.29470595 ,  5.44752274],[947.65535211 ,  4.7593997 ],[1021.51121852  ,  4.03057061],[1011.04537687 ,   4.02310175],[999.74521784 ,  3.99879389],[999.99900454  , 4.00005048],[1000.00047927   , 3.9999974 ],[999.99998606 ,  4.        ],[1000.00000017  ,  4.        ],[1000.  ,  4.]]
#constants = [[1.70716062,1.70716062,1.70716062],[2.05710384,2.64393153,1.70716062],[3.01684339,3.30562203,1.70716062],[8.64841707, 5.78859566, 1.70716062],[8.64995476, 5.77977766, 1.70716062],[8.6561055,  5.74450567, 1.70716062],[8.68070845, 5.60341771, 1.70716062],[8.75998067, 5.24648482, 1.70716062],[8.82088279, 5.04252095, 1.70716062],[8.9075418,  4.80105919, 1.70716062],[9.01861382, 4.52770821, 1.70716062],[9.01643462, 4.57224507, 1.70716062],[9.08256351, 4.67516206, 1.70716062],[9.3315203,  4.86929737, 1.70716062],[9.56764629, 4.97824891, 1.70716062],[9.75718046, 4.98552259, 1.70716062],[10.00479845,  4.78726893,  1.70716062],[10.08623312,  4.79331857,  1.70716062],[10.12336612,  4.8027204,   1.70716062],[10.12460006,  4.80365547,  1.70716062],[10.12456437,  4.80374017 , 1.70716062],[10.12455064,  4.80374721,  1.70716062],[10.12455064,  4.80374721,  1.70716062]]

#times = []

#print(len(constants))
#sys.exit()

#print(constants)
#for k_num in range(0,things):
#	#print('what?')
#	alter = pd.read_csv('./input_file.csv',header=None)
#	alter.iloc[33,1] = 'FALSE'
#	alter.iloc[34,1] = 'FALSE'
#	alter.iloc[35,1] = 'TRUE'
#	alter.iloc[40,1] = 'FALSE'
#	alter.iloc[17,1] = exp_data_folder+'iter_'+str(k_num)#'./analysis/CO_O2_analysis/CO_O2_fit_productobj_folder/fitting/iter_'+str(k_num)
#	#alter.iloc[17,1] = './CO_O2_fit_productobj_3rd_folder/fitting/iter_'+str(k_num)
#	value = 1
#	current_rate = 46
#		
#	for kz in range(0,1):
#		if value == 1:
#			alter.iloc[current_rate,value] = constants[k_num][kz]
#			value = 2
#		elif value == 2:
#			value = 2
#			if str(alter.iloc[current_rate,value]) != 'nan':
#				alter.iloc[current_rate,value] = constants[k_num][kz]
#			else:
#				pass
#			current_rate += 1
#			value = 1
#		
#	alter.to_csv('./input_file.csv',header=None,index=False)
#	#sys.exit()	
#	base_run()
#		#sys.exit()
constants = []
prop_mech = []
times = []
act_mech = []
generate_gif(reactants, syn_data_folder, exp_data_folder, num_in_folder, constants, prop_mech, times, act_mech)