
def vary_Input(variable_type=None,variableToChange=None, newValue=0, input_file='./input_file.csv'):
	df1 = pd.read_csv(input_file,header = None)
	a_vary_input = df1.where(df1==variableToChange).dropna(how='all').dropna(axis=1)
	cellCol = a_vary_input.columns[0]
	cellRow = a_vary_input.index[0]
	
	if variable_type == ('intensity' or 'surface') and variableToChange != None:
		df1.iloc[cellRow+1,cellCol] = newValue
	elif variable_type == 'delay' and variableToChange != None:
		df1.iloc[cellRow+2,cellCol] = newValue
	elif variable_type == 'Temperature':
		df1.iloc[cellRow,1] = newValue
	elif variableToChange == 'Output Folder Name':
		df1.iloc[cellRow,1] = newValue
	else:
		print('Variable not found in input file.')
		print('If not a user error, please contact developers.')
	df1.to_csv(input_file,header=None,index=False)