		if 'thermo equations' in reac_input.keys():
			if len(reac_input['thermo equations']) > 0:
				thermoConstraints = True

			else:
				thermoConstraints = False

		else:
			thermoConstraints = False


		if thermoConstraints == True:
			thermoReactions = {}
			thermoStoich = {}
			for znum,z in enumerate(reac_input['thermo equations']):
				thermoReactions[znum] = []
				thermoStoich[znum] = []
				tempValue = z.split('+')
				for j in tempValue:
					if j.find('*') != -1:
						thermoSplit = j.split('*')
						thermoSplit[1] = thermoSplit[1].replace('r','')
						thermoReactions[znum].append(int(thermoSplit[1]))
						thermoSplit[0] = thermoSplit[0].replace('(','')
						thermoSplit[0] = thermoSplit[0].replace(')','')
						thermoStoich[znum].append(float(thermoSplit[0]))
					else:
						j = j.replace('r','')
						thermoReactions[znum].append(int(j))
						thermoStoich[znum].append(1)


		if thermoConstraints == True:
			
			thermoReactants = {}
			thermoProducts = {}
			thermoReactantsValue = {}
			thermoProductsValue = {}

			for jnum,jval in enumerate(thermoReactions):
				thermoReactants[jnum] = []
				thermoProducts[jnum] = []
				thermoReactantsValue[jnum] = []
				thermoProductsValue[jnum] = []

				for znum,zval in enumerate(thermoReactions[jval]):

					if rev_irr[jval-1] != 1:
						print('')
						print('To set a thermodynamic constraints, all reactions must be reversible.')
						print('Currently, elementary reaction '+str(jval)+' is irreversible and included in the thermo. equation.')
						sys.exit()

					else:
						for zen in range(0,necessary_values['molecules_in_gas_phase']):
							print(rate_array[zval-1])
							if rate_array[zval-1][zen] == 1:

								thermoProductsValue[jnum].append(necessary_values['reactants'][zen])						
								
								if thermoStoich[jnum][znum] != 1:
									thermoProducts[jnum].append('('+str(thermoStoich[jnum][znum])+')*'+necessary_values['reactants'][zen])
									
								else:
									thermoProducts[jnum].append(necessary_values['reactants'][zen])	

							elif rate_array[zval-1][zen] == -1:

								thermoReactantsValue[jnum].append(necessary_values['reactants'][zen])
								if thermoStoich[jnum][znum] != 1:
									thermoReactants[jnum].append('('+str(thermoStoich[jnum][znum])+')*'+necessary_values['reactants'][zen])
								else:
									thermoReactants[jnum].append(necessary_values['reactants'][zen])		
				print(thermoReactantsValue)
				print(thermoProductsValue)
				print(thermoStoich)

				overallReactionString = ''
			
				for znum in range(0,len(thermoReactants[jnum])):
					overallReactionString = overallReactionString+thermoReactants[jnum][znum]
					if znum != len(thermoReactants[jnum])-1:
						overallReactionString = overallReactionString + ' + '

				overallReactionString = overallReactionString + ' <-> '
				
				for znum in range(0,len(thermoProducts[jnum])):	
					overallReactionString = overallReactionString + thermoProducts[jnum][znum]
					if znum != len(thermoProducts[jnum])-1:
						overallReactionString = overallReactionString + ' + '
				print('Thermodynamic Constraint Equation:')
				print(overallReactionString)
				print('')

				if reactor_kinetics_input['thermo values'][jnum] == 'none':
					freeEnergyValue = 0

					try:
						for znum in range(0,len(thermoReactantsValue[jnum])):
							calcValue = (-1)*(thermoStoich[jnum][znum])*molecularProperties(thermoReactantsValue[jnum][znum],'freeEnergy',temperature=reac_input['Reactor Temperature'])
							freeEnergyValue += calcValue

						for znum in range(0,len(thermoProductsValue[jnum])):

							calcValue = (thermoStoich[jnum][znum])*molecularProperties(thermoProductsValue[jnum][znum],'freeEnergy',temperature=reac_input['Reactor Temperature'])
							
							freeEnergyValue += calcValue
					except:
						print('Failed to calculate free energy of reaction automatically. Please enter free energy manually.')
						sys.exit()

					reactor_kinetics_input['thermo values'][jnum] = freeEnergyValue
					print('Free energy calculated automatically: ')
					print(reactor_kinetics_input['thermo values'][jnum])
					print('')
					

				else:
					print('Free energy entered manually: ')
					print(reactor_kinetics_input['thermo values'][jnum])
					print('')
