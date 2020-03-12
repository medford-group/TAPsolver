import pandas as pd 
import numpy as np
import sys
import os

names = ['Inert-1','O216','O16O18','O218','Inert-2']

for j in range(1,6):
	
	data = pd.read_csv(str(j)+'.csv')
	print(data.shape[1])

	for k in range(1,data.shape[1]):
		try:
			os.mkdir('./flux_data/pulse'+str(k))
		except:
			pass
		newFile = pd.concat([data.iloc[:,0],data.iloc[:,k]],axis = 1)

		newFile.to_csv('./flux_data/pulse'+str(k)+'/'+names[j-1]+'.csv', index=False)
		print(names[j-1])
	#to_csv

	#print(data)

	#print(type(data))