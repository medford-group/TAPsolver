def bulk_structure(reactor_setup,iterDict,baseName=None):

		# generate csv submission calculation

		directoryGenerator('./'+baseName)
		copyfile(reactor_setup,'./'+baseName+'/'+reactor_setup)

		if baseName == None:
			print('No variables are changed. No generation performed.')
			sys.exit()

		#headers = []
		headers = ['directory']

		generalVariables = list(iterDict.keys())

		for j in iterDict.keys():
			if j == 'mechanisms':
				headers.append(j)
			else:
				for k in iterDict[j].keys():
					headers.append(k+'_'+j)

		def iterativeGeneration(newNpArray,directoryName,currentList,currentkey,currentSubKey):

			if list(iterDict.keys())[currentkey] == "mechanisms":

				for jnum,j in enumerate(list(iterDict[list(iterDict.keys())[currentkey]].keys())):

					newdirectoryName = directoryName+str(jnum)
					newList = [j]

					copyfile(iterDict["mechanisms"][j],'./'+baseName+'/'+iterDict["mechanisms"][j])
					# Construct new base input_file
					input_construction(reactor_name='./'+reactor_setup,reaction_name=iterDict["mechanisms"][j],new_name='./input_file.csv')
					newNpArray = iterativeGeneration(newNpArray,newdirectoryName,newList,currentkey+1,0)

			else:
				#tempInput = pd.read_csv('./input_file.csv',header=None)
				for jnum,j in enumerate(iterDict[list(iterDict.keys())[currentkey]][list(iterDict[list(iterDict.keys())[currentkey]].keys())[currentSubKey]]):
					newdirectoryName = directoryName+str(jnum)				
					newdirectoryName = directoryName+str(jnum)
					newList = currentList[:]
					newList.append(j)

					tempInput = pd.read_csv('./input_file.csv',header=None)
					currentSpecies = list(iterDict[list(iterDict.keys())[currentkey]])[currentSubKey]
					print(tempInput)
					if list(iterDict.keys())[currentkey] == "surface":
						currentSurface = list(tempInput.iloc[21,1:])
						tempInput.iloc[22,1+currentSurface.index(currentSpecies)] = float(j)
						tempInput.to_csv('./input_file.csv',header=None,index=False)
					elif list(iterDict.keys())[currentkey] == "intensity":
						currentGasses = list(tempInput.iloc[16,1:])
						tempInput.iloc[17,1+currentGasses.index(currentSpecies)] = float(j)
						tempInput.to_csv('./input_file.csv',header=None,index=False)
					elif list(iterDict.keys())[currentkey] == "time":
						currentGasses = list(tempInput.iloc[16,1:])
						tempInput.iloc[18,1+currentGasses.index(currentSpecies)] = float(j)
						tempInput.to_csv('./input_file.csv',header=None,index=False)
					if currentSubKey+1 != len(list(iterDict[list(iterDict.keys())[currentkey]].keys())):
						newNpArray = iterativeGeneration(newNpArray,newdirectoryName,newList,currentkey,currentSubKey+1)
					elif currentkey+1 != len(list(iterDict.keys())):
						newNpArray = iterativeGeneration(newNpArray,newdirectoryName,newList,currentkey+1,0)
					else:
						directoryNameList = [newdirectoryName]
						newList = directoryNameList+newList
						newNpArray = np.vstack((newNpArray,newList))
						directoryGenerator(newdirectoryName)
						tempInput.to_csv('./input_file.csv',header=None,index=False)
						
						copyfile('./running_simulation.py',newdirectoryName+'/running_simulation.py')
						copyfile('./run.sh',newdirectoryName+'/run.sh')
						copyfile('./input_file.csv',newdirectoryName+'/input_file.csv')

			return newNpArray
						
		dataStructure = np.asarray(headers)
		dataValues = np.asarray(headers[2:])

		dataStructurenew = iterativeGeneration(dataStructure,'./'+baseName+'/',[],0,0)
		df = pd.DataFrame(dataStructurenew)
		df.to_csv('./'+baseName+'/'+'directoryList.csv',header=None,index=False)

		valueList = []

		for j in list(iterDict.keys())[1:]:
			for k in iterDict[j].keys():
				valueList.append(iterDict[j][k])

		dataValues = np.vstack((list(dataValues),valueList))

		df2 = pd.DataFrame(dataValues)

		df2.to_csv('./'+baseName+'/'+'changedVariables.csv',header=None,index=False)
