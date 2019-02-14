import pandas as pd
import sys
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import re
from reaction_list_odes import deriv_and_const


def R_RRM_func(input_reactants,current_path,data_folder):

	path = current_path+'/'+data_folder+'/RRM_results/'
	try:  
		os.mkdir(path)
	except OSError:  
		print ("Creation of the directory %s failed" % path)
		pass
	else:  
		print ("Successfully created the directory %s " % path)
	
	#for j in r_const:
	#	r_const[j] = Constant(r_const[j])

	#user_data = pd.read_csv('./input_file.csv',header=None)
	
	#reactor_info = user_data.iloc[2:18,:] 
	#feed_surf_info = user_data.iloc[21:27,:]
	#data_storage = user_data.iloc[29:35,:]
	#reaction_info = user_data.iloc[38:,:]
	

	user_data = pd.read_csv('./input_file.csv',header=None)

	rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
	rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
	rows_3, cols_3 = np.where(user_data == 'Data_Storage_Options')
	rows_4, cols_4 = np.where(user_data == 'Reaction_Information')

	reactor_info = user_data.iloc[(2+rows_1[0]):(rows_2[0]-1),:] 
	feed_surf_info = user_data.iloc[2+rows_2[0]:rows_3[0]-1,:]
	data_storage = user_data.iloc[2+rows_3[0]:rows_4[0]-1,:]
	reaction_info = user_data.iloc[2+rows_4[0]:,:]

	reactor_kinetics_input = {}

	for k in range(0,len(feed_surf_info.index)):
		try:
			reactor_kinetics_input[feed_surf_info.iloc[k,5]] = float(feed_surf_info.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[feed_surf_info.iloc[k,5]] = feed_surf_info.iloc[k,1]



	reactor_info = reactor_info.T
	a = reactor_info.iloc[5]
	reactor_info.iloc[0] = a
	reactor_info = reactor_info.iloc[:2,:]
	
	feed_surf_info = feed_surf_info.T
	a = feed_surf_info.iloc[5]
	feed_surf_info.iloc[0] = a
	feed_surf_info = feed_surf_info.iloc[:2,:]
	
	reactor_info = pd.concat([reactor_info,feed_surf_info], axis=1)
	
	data_storage = data_storage.T
	a = data_storage.iloc[5]
	data_storage.iloc[0] = a
	data_storage = data_storage.iloc[:2,:]
	
	reactor_info = pd.concat([reactor_info,data_storage], axis=1)
	
	reaction_info = reaction_info.T
	a = reaction_info.iloc[5]
	reaction_info.iloc[0] = a
	reaction_info = reaction_info.iloc[:2,:]
	
	reactor_info = pd.concat([reactor_info,reaction_info], axis=1)
	
	reactor_info.to_csv(data_folder+'/RRM_params.csv',header=None)
	
	#print(reactor_info)
	
	command = 'Rscript'
	
	#path2script = './simpleRRM.R'#'./yongeRRM.R'
	path2script = './yongeRRM.R'#'./yongeRRM.R
	
	feed = 1
	reactants = input_reactants
	in_reactants = reactants
	print(reactants)
	reactants = str(reactants)
	reactants = reactants.replace('[','')
	reactants = reactants.replace(']','')
	reactants = reactants.replace("'","")
	reactants = reactants.replace(' ','')
	
	#sys.exit()
	
	arg1 = current_path+'/'+data_folder+'/'#"/home/adam/research_medford/python_code/tap_code/csv_input_sim/eley_eluc_folder/"
	arg2 = reactants
	arg3 = str(feed)
	
	pass_arg = []
	pass_arg.append(arg1)
	pass_arg.append(arg2)
	pass_arg.append(arg3)
	#sys.exit()
	cmd = [command, path2script,arg1,arg2,arg3] # ... <- the arguments that you need to call to get the script to run properly
	
	subprocess.run(cmd) #Will just save the output data as a file
	

	textstr2 = 'RRM Model:'
	
	for k in range(len(in_reactants)):
		new_model = ''
		new_model += 'R_'+in_reactants[k]+' = '
		new_model += '\u03A8'+'_'+str(0)+' + '
		for j in range(len(in_reactants)):
			new_model += '\u03A8'+'_'+str(j+1)+'*C_'+''+in_reactants[j]+' + '
	
		for j in range(len(in_reactants)):
			new_model += '\u03A8'+'_'+str(len(in_reactants)+j+1)+'*U_'+''+in_reactants[j]+' + '
	
		for j in range(len(in_reactants)):
			for k in range(len(in_reactants)):
				new_model += '\u03A8'+'_'+str((j+2)*len(in_reactants)+k+1)+'*C_'+''+in_reactants[j]+'*U_'+''+in_reactants[k]+' + '
	
		for j in range(len(in_reactants)):
			for k in range(j,len(in_reactants)):
				#print((j-1)*len(in_reactants)+2*len(in_reactants)+len(in_reactants)*len(in_reactants)+k+1)
				new_model += '\u03A8'+'_'+str((j)*len(in_reactants)+2*len(in_reactants)+len(in_reactants)*len(in_reactants)+k+1-j)+'*U_'+''+in_reactants[j]+'*U_'+''+in_reactants[k]
				if j == (len(in_reactants)-1) and k == (len(in_reactants)-1):
					pass
				else:
					new_model += ' + '
		#print(new_model)
		#print()
		
		textstr2 = '\n'.join((textstr2,new_model))
	
	#	print()
	#	print()
	

	#total_u = 0
	#for j_k in range(0,len(in_reactants)):
	#	total_u += len(in_reactants)-j_k
	#total_u += len(in_reactants)*len(in_reactants)+len(in_reactants)*len(in_reactants)
	#print(total_u)
	#sys.exit()
	
	df_time = pd.read_csv(current_path+'/'+data_folder+'/flux_data/Inert.csv',header=None).iloc[:,0]

	rows = in_reactants.copy()
	for k in in_reactants.copy():
		rows.append('U_'+k)

	RRM_data = {}

	for k in in_reactants:
		RRM_data[k] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+k+'_reactivities.csv',header=None)

	y_proc_data = {}
	for k in range(int(reactor_kinetics_input['number_of_pulses'])):
		y_proc_data[k+1] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+str(k+1)+'_y_proc.csv',header=None).iloc[1:,:]
	#new_frame = y_proc_data[0+1].iloc[:,1+0+len(in_reactants)].copy()
	#print(new_frame)
	#sys.exit()
	
	#fig2, ax2 = plt.subplots()
	#for k in range(1,len(y_proc_data[1].columns)-3):
	#	#pd.to_numeric(new_frame,errors='coerce')
	#	ax2.plot(df_time,pd.to_numeric(y_proc_data[1].iloc[:,k]))
	#plt.show()

	#step through the pulses
	for z in range(int(reactor_kinetics_input['number_of_pulses'])):
		f, axarr = plt.subplots(len(in_reactants)*2, len(in_reactants))
		f = plt.gcf()
		f.set_size_inches(3*(len(in_reactants)), 1.5*len(rows))

		try:
			for ax, col in zip(axarr[0], in_reactants):
				#new_in = r'$R_%d$' % col
				ax.set_title(r'$\partial$'+'('+col+')'+' /')
				for ax, row in zip(axarr[:,0], rows):
					ax.set_ylabel('/ '+r'$\partial$'+'('+row+')', rotation=90, size='large')
		except TypeError:
			axarr[0].set_title(cols[0])
			for ax, row in zip(axarr[:], rows):
				ax.set_ylabel('/ '+r'$\partial$'+row, rotation=90, size='large')


		#step through the numerator
		for k in range(len(in_reactants)):
			#step through the gas denominator
			for j in range(len(in_reactants)):
				
				try:
					test = axarr[j,k]
				except IndexError:
					test = axarr[j]

				df_merged = pd.DataFrame(index=y_proc_data[1].index,columns=range(1))#.shape[0]

				df_merged.iloc[:] = float(0)#float(RRM_data[in_reactants[k]].iloc[z+1,j+2])
				
				#RRM_data[k] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+k+'_reactivities.csv',header=None)

				#print('dR_'+in_reactants[k]+' / dC_'+in_reactants[j],end=' = ')
				#print('\u03A8'+'_'+str((j+1)),end=' + ')

				for j_nu in range(0,len(in_reactants)):
					new_frame = y_proc_data[z+1].iloc[:,1+j_nu]#.copy()
					
					new_frame = pd.to_numeric(new_frame,errors='coerce')
					
					new_frame *= float(RRM_data[in_reactants[k]].iloc[z+1,2*len(in_reactants)+(k)*len(in_reactants)+2+j_nu])

					df_new = new_frame.to_frame()
					df_new.columns = [0]
					df_merged = df_merged.add(df_new,fill_value=0)
					#df_merged = df_merged + new_frame.to_frame()
					#df_merged[0] = df_merged.iloc[0] + new_frame.to_frame().iloc[0]
					#print(df_merged)
					#sys.exit()
					#df_merged.iloc[0] = df_merged.iloc[0] + new_frame#+ RRM_data[in_reactants[k]].iloc[z+1,2*len(in_reactants)+(k)*len(in_reactants)+2+j_nu]*y_proc_data[z+1].iloc[:,j_nu]

				#sys.exit()
				#print(df_merged)
				#sys.exit()

				#for j_nu in range(0,len(in_reactants)):
				#	print('\u03A8'+'_'+str(2*len(in_reactants)+(k)*len(in_reactants)+1+j_nu)+' * U_'+in_reactants[j_nu],end=' + ')
				test.plot(df_time, df_merged[0])

			#step through the surface denominator
			for j in range(len(in_reactants)):
				try:
					test = axarr[j+len(in_reactants),k]
				except IndexError:
					test = axarr[j+len(in_reactants)]

				df_merged = pd.DataFrame(index=y_proc_data[1].index,columns=range(1))#.shape[0]

				df_merged.iloc[:] = 0#float(RRM_data[in_reactants[k]].iloc[z+1,len(in_reactants)+j+2])



				#aditional terms
				#for step_1 in range(0,len(in_reactants)):

				print()
				# Graphs for dU portion
				#for j_nu in range(0,len(in_reactants)):
				curr_value = k + 2*len(in_reactants) 
				print(curr_value)
				for step_2 in range(0,k):
					new_frame = y_proc_data[z+1].iloc[:,1+step_2+len(in_reactants)].copy()
					#print(new_frame)
					
					new_frame = pd.to_numeric(new_frame,errors='coerce')

					new_frame *= float(RRM_data[in_reactants[k]].iloc[z+1,curr_value+1])
					#print(float(RRM_data[in_reactants[k]].iloc[z+1,curr_value+1]))

					df_new = new_frame.to_frame()
					df_new.columns = [0]
					df_merged = df_merged.add(df_new,fill_value=0)

					curr_value += (len(in_reactants)-step_2-1)

					df_merged.iloc[0] = df_merged.iloc[0] + new_frame
					print(curr_value)
				#if k > 1:
				#	sys.exit()
				curr_value = k + 2*len(in_reactants) + len(in_reactants)*len(in_reactants) + 1 #step_1
				for step_3 in range(k,len(in_reactants)):#curr_value+
					new_frame = y_proc_data[z+1].iloc[:,1+step_3].copy()
					new_frame = pd.to_numeric(new_frame,errors='coerce')
					#new_frame *= float(RRM_data[in_reactants[k]].iloc[z+1,curr_value+1])
					curr_value += 1
					df_new = new_frame.to_frame()
					df_new.columns = [0]
					df_merged = df_merged.add(df_new,fill_value=0)
					#print('\u03A8'+'_'+str(curr_value)+' * U_'+in_reactants[step_3],end=' + ')
					df_merged.iloc[0] = df_merged.iloc[0] + new_frame
				#sys.exit()
				#print(df_time)
				#print(df_merged[0])
				#print(type(df_time))
				#print(type(df_merged[0]))
				#sys.exit()
				test.plot(df_time, df_merged)
				
				
				####Generate the derivative equations
				#print('dR_'+in_reactants[k]+'/dU_'+in_reactants[j],end=' = ')
				#print('\u03A8'+'_'+str(len(in_reactants)+(j+1)),end=' + ')
				#curr_value = k + 2*len(in_reactants) + len(in_reactants)*len(in_reactants) + 1 #step_1
				#for j_nu in range(0,len(in_reactants)):
				#	print('\u03A8'+'_'+str(2*len(in_reactants)+(j_nu)*len(in_reactants)+1+j)+' * C_'+in_reactants[j_nu],end=' + ')


				#for step_2 in range(0,k):
				#	print('\u03A8'+'_'+str(curr_value)+' * U_'+in_reactants[step_2],end=' + ')
					
			
				#for step_3 in range(k,len(in_reactants)):#curr_value+
				#	print('\u03A8'+'_'+str(curr_value)+' * U_'+in_reactants[step_3],end=' + ')
				#	curr_value += 1
				#print()
				#print()


		props = dict(boxstyle='round', facecolor='white', alpha=0.4)
		plt.gcf().text(0.5, 0.05, textstr2, fontsize=8, bbox=props,ha='center')
		f.tight_layout()	
		f.subplots_adjust(bottom=.2,left=0.3,top=0.9,right=0.7)
		f.suptitle("RRM Jacobian", fontsize=14)
		#f.text(0.15, 0.6, 'Derivative Denominator', va='center', rotation='vertical', fontsize=14)
		
		plt.savefig(current_path+'/'+data_folder+'/graphs/RRM_p'+str(z+1)+'.png')#current_path+'/'+data_folder+
		plt.show()

def petal_plots_RRM(folder_location):

	print('good pass')