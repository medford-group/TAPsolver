def generateFolder(path_name):
	try:  
		os.mkdir(path_name)
	except OSError:  
		pass
	else:  
		pass
	
def storeSens(yes_no,output_file,gasses,legend_ref):	
	if yes_no == True:
		try:
			os.makedirs('./'+output_file+'_folder/sensitivity_'+output_file)
		except OSError:
			pass
		for k_sens_folder in range(gasses):	
			try:
				os.makedirs('./'+output_file+'_folder/sensitivity_'+output_file+'/'+legend_ref[k_sens_folder])
			except OSError:
				pass

def storeDataFunc(yes_no,output_file):
	if yes_no == True:
		try:
			os.makedirs('./'+output_file+'_folder')
		except OSError:
			pass

def directory_generator(path_name):
	try:  
		os.mkdir(path_name)
	except OSError:  
		pass
	else:  
		pass