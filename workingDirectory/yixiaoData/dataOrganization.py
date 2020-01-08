import pandas as pd 
import numpy as np 
import sys
import os

weights = [0,1,2,3,4]
for j in weights:
	weight = j
	
	speciesList = ['O216','O218','O16O18','Ar','He']
	
	# All file options
	Mn = 'Mn_t=2s.xlsx'
	
	# Define specific spreadsheet
	xls = pd.ExcelFile(str(weight)+Mn)
	
	# Read specific spreadsheets
	O216 = pd.read_excel(xls, '2')
	O218 = pd.read_excel(xls, '4')
	O16O18 = pd.read_excel(xls, '3')
	Ar = pd.read_excel(xls, '1')
	He = pd.read_excel(xls, '5')
	
	# Make folder (test existence first)
	if not os.path.exists('./'+str(weight)+'_folder'):
		os.mkdir('./'+str(weight)+'_folder')
	
	if not os.path.exists('./'+str(weight)+'_folder/flux_data'):
		os.mkdir('./'+str(weight)+'_folder/flux_data')
	
	#for k_str in speciesList:
	#	if not os.path.exists('./'+str(weight)+'_folder/flux_data/'+k_str):
	#		os.mkdir('./'+str(weight)+'_folder/flux_data/'+k_str)

	np.savetxt('./'+str(weight)+'_folder/flux_data/'+'Inert-1.csv',Ar.iloc[1:5000,2:],delimiter=",")
	np.savetxt('./'+str(weight)+'_folder/flux_data/'+'Inert-2.csv',He.iloc[1:5000,2:],delimiter=",")
	np.savetxt('./'+str(weight)+'_folder/flux_data/'+'O216.csv',O216.iloc[1:5000,2:],delimiter=",")
	np.savetxt('./'+str(weight)+'_folder/flux_data/'+'O218.csv',O218.iloc[1:5000,2:],delimiter=",")
	np.savetxt('./'+str(weight)+'_folder/flux_data/'+'O16O18.csv',O16O18.iloc[1:5000,2:],delimiter=",")