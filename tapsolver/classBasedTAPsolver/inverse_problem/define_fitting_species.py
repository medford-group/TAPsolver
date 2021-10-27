
		if reactor_kinetics_input['Fit Parameters'].lower() == 'true' or (reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' and (sens_type == 'total' or sens_type == 'initial')) or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
			speciesString = ''
			for k in range(0,int(necessary_values['molecules_in_gas_phase'])):
				if k != int(necessary_values['molecules_in_gas_phase']):
					speciesString+='1,'
				elif int(reac_input['Number of Inerts']) != 0:
					speciesString+='1'
			for k in range(0,int(reac_input['Number of Inerts'])):
				if k != int(reac_input['Number of Inerts']):
					speciesString+='0,'
				else:
					speciesString+='0'

		elif reactor_kinetics_input['Fit Inert'].lower() == 'true':
			speciesString = ''
			for k in range(0,int(necessary_values['molecules_in_gas_phase'])):
				if k != int(necessary_values['molecules_in_gas_phase']):
					speciesString+='0,'
				elif int(reac_input['Number of Inerts']) != 0:
					speciesString+='0'
			for k in range(0,int(reac_input['Number of Inerts'])):
				if k != int(reac_input['Number of Inerts']):
					speciesString+='1,'
				else:
					speciesString+='1'

		else:
			speciesString = ''
			for k in range(0,int(necessary_values['molecules_in_gas_phase'])):
				if k != int(necessary_values['molecules_in_gas_phase']):
					speciesString+='0,'
				elif int(reac_input['Number of Inerts']) != 0:
					speciesString+='0'
			for k in range(0,int(reac_input['Number of Inerts'])):
				if k != int(reac_input['Number of Inerts']):
					speciesString+='0,'
				else:
					speciesString+='0'

		reactor_kinetics_input['Objective Species'] = speciesString
