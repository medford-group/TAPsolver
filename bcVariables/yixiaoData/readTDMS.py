from nptdms import TdmsFile
import numpy as np
import sys

#tdms_file = TdmsFile('./4MnNa2WO4/4Mn5Na2WO4_t=2s.tdms')
tdms_file = TdmsFile('./WO4/5WO4_t=2s.tdms')


#tdms_groups = tdms_file.groups()

#print(tdms_groups)

#tdms_Variables_1 = tdms_file.group_channels('Meta Data')

#print(tdms_Variables_1)

#MessageData_channel_1 = tdms_file.object('Meta Data','Item')

#print(MessageData_channel_1.data)

j = 0
k = 0
for group in tdms_file.groups():
	start = True
	print(group)
	j+=1

	for channel in tdms_file.group_channels(group):

		properties = channel.properties
		data = channel.data
		currentDataType = type(data)

		print(group)
		print(data)
			
		k+=1
		
		#print(group)
		#print(type(channel))

		
		if currentDataType == np.ndarray:
			print(data[::10])		
			if len(data) > 49999:
	
				if start == True:
					newTotalArray = data[::10]
					start = False
				
				else:
					try:
						newTotalArray = np.vstack((newTotalArray,data[::10]))
					except:
						newTotalArray = np.vstack((newTotalArray,data[::10][:-1]))
	try:
		np.savetxt('./'+group+'.csv', np.transpose(newTotalArray), delimiter=",")
		
	except:
		
		pass
