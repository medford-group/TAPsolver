import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys
from reaction_list_odes import deriv_and_const

def MKM_graphs(kin_params,reactions_in,fold_loc,graph_display):
	#print(kin_params)
	ke0 = 1
	gasses = 2
	
	folder_location = './'+fold_loc+'/'#./new_test_folder/'
	
	reactions = reactions_in#['CO + * <-> CO*','CO + O* <-> CO2','CO* + O* <-> CO2']
	reactions_input = reactions.copy()
	deriv_dict, reactants, rate_equations = deriv_and_const(reactions,gasses)
	
	
	deriv = len(reactants)

	textstr = 'Elementary Reactions:'
	for k_elem in reactions_input:
		textstr = '\n'.join((textstr,k_elem))
	
	textstr2 = 'Rate Expressions:'
	for k_reac in rate_equations:
		textstr2 = '\n'.join((textstr2,k_reac))
	
	cols = reactants[:gasses]
	
	rows = reactants#['CO','CO2','CO*','O*','*']
	
	df_time = pd.read_csv(folder_location+'*'+'.csv',header=None).iloc[:,0]
	
	#example_array = [['-1*(ke0*[*]*[CO*]) + -1*(ke1*[O*])','1*(ke1*[O*])'],['-1*(kd1)','0'],['0','1*(ke2*[O*])'],['-1*(ke1*[CO])','1*(ke1*[CO]) + 1*(ke2*[CO*])'],['-1*(ke0*[CO])','0']]
	
	f, axarr = plt.subplots(deriv, gasses)
	
	f = plt.gcf()
	f.set_size_inches(3*len(cols), 1.5*len(rows))
	#fig.savefig('test2png.png', dpi=100)
	
	
	try:
		for ax, col in zip(axarr[0], cols):
		    ax.set_title(col)
		    for ax, row in zip(axarr[:,0], rows):
	    		ax.set_ylabel(row, rotation=90, size='large')
	except TypeError:
		axarr[0].set_title(cols[0])
		for ax, row in zip(axarr[:], rows):
			ax.set_ylabel(row, rotation=90, size='large')
	
	
	for k_der in range(0,deriv):
		for j_spe in range(0,gasses):
			try:
				test = axarr[k_der,j_spe]
			except IndexError:
				test = axarr[k_der]
			test.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			try:
				new = deriv_dict[(j_spe,k_der)].split('+')#example_array[k_der][j_spe].split('+')
			except NameError:
				new = deriv_dict[(j_spe,k_der)]#example_array[k_der][j_spe]
			df_merged = pd.DataFrame(index=df_time.index,columns=range(1))
			df_merged.iloc[:] = 0
			for z in new:
				if '(' in z:
					result = (z.split('('))[1].split(')')[0]
					#result = re.search('%s(.*)%s' % ('(', ')'), z)
					result = result.split('*[')
					for jk_num,jk in enumerate(result.copy()):
						if "]" in jk:
							result[jk_num] = jk.replace("]",'')
					#print(result)
					if len(result) > 1:
						new_df = pd.read_csv(folder_location+'/'+result[1]+'.csv',header=None).iloc[:,1]
						if len(result) > 2:
							for k_cool in range(2,len(result)):
								#print(result[k_cool])
								df_new = pd.read_csv(folder_location+'/'+result[k_cool]+'.csv',header=None).iloc[:,1]
								new_df = new_df*df_new
						if '-' in z:
							new_df *= -kin_params[result[0].replace('k', 'K')]#ke0#str
						else:
							new_df *= kin_params[result[0].replace('k', 'K')]#kin_params[result[0]]#ke0
					elif (z.split('('))[1].split(')')[0] != '':
						new_df = pd.DataFrame(index=df_time.index,columns=range(1))
						if '-' in z:
							new_df.iloc[:] = -kin_params[result[0].replace('k', 'K')]#kin_params[result[0]]
						else:
							new_df.iloc[:] = kin_params[result[0].replace('k', 'K')]#kin_params[result[0]]#ke0
					else:
						new_df.iloc[:] = 0
				temp = df_merged 
				df_merged[0] = temp[0] + new_df
				
				#print(df_merged)
					#test.plot(df_time, df_merged)
			#print(df_merged)
			#sys.exit()
			props = dict(boxstyle='round', facecolor='white', alpha=0.0)
			#axarr[k_der, j_spe].set_xlabel(example_array[k_der][j_spe])
			if deriv_dict[(j_spe,k_der)] != '0':#example_array[k_der][j_spe] != '0':
				result = (deriv_dict[(j_spe,k_der)].split('('))[1].split(')')[0]
				if result != '':
					pass
					#if '+' in deriv_dict[(j_spe,k_der)]:#example_array[k_der][j_spe]:
					#	new_label = deriv_dict[(j_spe,k_der)].split('+')#example_array[k_der][j_spe].split('+')
					#	#test.text(0.5, 0.77, new_label[0],str('+')+new_label[1], transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')
					#	test.text(0.5, -0.3, '\n'.join((new_label[0],str('+')+new_label[1])), transform=test.transAxes, fontsize=7,verticalalignment='bottom', bbox=props,ha='center')
					#else:	
					#	test.text(0.5, -0.3, deriv_dict[(j_spe,k_der)], transform=test.transAxes, fontsize=7,verticalalignment='bottom', bbox=props,ha='center')#, bbox=props
					#	#test.text(0.5, 0.77, '\n'.join(('Deriv:',deriv_dict[(j_spe,k_der)])), transform=test.transAxes, fontsize=8,verticalalignment='top', bbox=props,ha='center')
			else:
				df_merged = pd.DataFrame(index=df_time.index,columns=range(1))
				df_merged.iloc[:] = 0
			
			test.plot(df_time, df_merged)

	props = dict(boxstyle='round', facecolor='white', alpha=0.4)
	plt.gcf().text(0.75, 0.85, textstr, fontsize=8, bbox=props,ha='left')
	plt.gcf().text(0.5, 0.05, textstr2, fontsize=8, bbox=props,ha='center')
	f.tight_layout()	
	f.subplots_adjust(bottom=.2,left=0.3,top=0.9,right=0.7)
	f.suptitle("MKM Jacobian", fontsize=14)
	plt.savefig(fold_loc+'/../graphs/Jacob_MKM.png')
	
	if graph_display == 'true':
		plt.show()
	else:
		plt.clf()
		plt.close()
	
				#try:
				#	result = re.search('[(.*)]',z)
				#	print(result.group(0))
				#	print('Test')
				#except TypeError:
				#	print('Failed')
				#	pass
	
				#if '-' in new:
				#	print('negative')
				#else:
				#	print('positive')

def petal_plots_exp(folder_location):

	print('good pass')