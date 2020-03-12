import pandas as pd 
import sys

species_list = ['Inert-1','Inert-2','O16O18','O216','O218']

for j in species_list:
	dataFrame1 = pd.read_csv('./flux_data/pulse1/'+j+'.csv')
	dataFrame2 = dataFrame1.copy()
	
	
	data1 = dataFrame1.iloc[:899,0]
	data1 = data1.reset_index(drop=True)
	
	data2 = dataFrame2.iloc[2099:2999,1]
	data2 = data2.reset_index(drop=True)
	
	test = pd.concat([data1,data2],axis=1,ignore_index=True)

	test.to_csv('./'+j+'.csv',index=False,header=False)