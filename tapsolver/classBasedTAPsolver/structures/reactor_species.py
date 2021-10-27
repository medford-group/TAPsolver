
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved


import pandas as pd
import numpy as np
import math as mp
from reactor_species import define_gas, define_adspecies

class reactor_species():

	"""
	
	This class acts as a container for all of the assumed initial conditions of an experiment.
	
	Args:
		
		gasses (dict of ): The  

	"""

	def __init__(self, inert_diffusion = 16, catalyst_diffusion = 16, reference_temperature = 385.6, reference_mass = 40, temerature = 385.6):
		self.gasses = {}
		self.inert_gasses = {}
		self.adspecies = {}
		self.inert_diffusion = inert_diffusion
		self.catalyst_diffusion = catalyst_diffusion
		self.reference_temperature = reference_temperature
		self.reference_mass = reference_mass
		self.temperature = temerature
		self.advection = 0
		self.reference_pulse_size = 1
		
	def add_gas(self,name='', define_gas_data = define_gas):
		def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
			return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))

		if name not in self.gasses:
			self.gasses[name] = define_gas_data
			self.gasses[name].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.gasses[name].mass)
			self.gasses[name].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.gasses[name].mass)
		else:
			print('Gas already defined in dictionary.')

	def add_inert_gas(self,name='', define_gas_data = define_gas):
		def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
			return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))

		if name not in self.gasses:
			self.inert_gasses[name] = define_gas_data
			self.inert_gasses[name].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.inert_gasses[name].mass)
			self.inert_gasses[name].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.inert_gasses[name].mass)
		else:
			print('Gas already defined in dictionary.')

	def add_adspecies(self,name='', define_adspecies_data = define_adspecies):
		
		if name not in self.adspecies:
			self.adspecies[name] = define_adspecies_data
		else:
			print('Gas already defined in dictionary.')

	@property
	def inert_diffusion(self):
		return self._inert_diffusion
	@inert_diffusion.setter
	def inert_diffusion(self,value):		
		def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
			return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	

		if value <= 0:
			raise ValueError("Reference inert diffusion must be positive (non-negative")
		self._inert_diffusion = value
		for j in self.gasses:
			self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.gasses[j].mass)
		for j in self.inert_gasses:
			self.inert_gasses[j].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self._inert_gasses[j].mass)

	@property
	def catalyst_diffusion(self):
		return self._catalyst_diffusion
	@catalyst_diffusion.setter
	def catalyst_diffusion(self,value):		
		def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
			return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	

		if value <= 0:
			raise ValueError("Reference catalyst diffusion must be positive (non-negative")
		self._catalyst_diffusion = value
		for j in self.gasses:
			self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.gasses[j].mass)
		for j in self.inert_gasses:
			self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self._inert_gasses[j].mass)

	@property
	def reference_temperature(self):
		return self._reference_temperature
	@reference_temperature.setter
	def reference_temperature(self,value):		
		def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
			return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	

		#if type(value) != Constant: 
		if value <= 0:
			raise ValueError("Temperature must be positive (non-negative")
		self._reference_temperature = value
		for j in self.gasses:
			self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self._reference_temperature, self.temperature, self.gasses[j].mass)
			self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self._reference_temperature, self.temperature, self.gasses[j].mass)

	@property
	def reference_mass(self):
		return self._reference_mass
	@reference_mass.setter
	def reference_mass(self,value):		
		def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
			return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	
		if value <= 0:
			raise ValueError("Species mass must be positive (non-negative")
		self._reference_mass = value
		for j in self.gasses:
			self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self._reference_mass, self.reference_temperature, self.temperature, self.gasses[j].mass)
			self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self._reference_mass, self.reference_temperature, self.temperature, self.gasses[j].mass)

	@property
	def temperature(self):
		return self._temperature
	@temperature.setter
	def temperature(self,value):
		if value == (float or int):		
			def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
				return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	

			if value <= 0:
				raise ValueError("Temperature must be positive (non-negative")
			self._temperature = value
			for j in self.gasses:
				self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self.reference_temperature, self._temperature, self.gasses[j].mass)
				self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self.reference_temperature, self._temperature, self.gasses[j].mass)
		else:
			self._temperature = value
	## Total diffusion matrix 

	## initializeVariableDictionaries - > in forward problem?

	## Scaling outlet flux - > in forward problem?

	## Adjusting inert gas diffusions -> in forward problem?

	## Pulse functions  - > in forward problem?

	## Surface functions  - > in forward problem?

	## 
	

	#@property
	#def inert_diffusion(self):
	#	#if self._inert_diffusion
	#	return self._inert_diffusion
	#@inert_diffusion.setter
	#def inert_diffusion(self,value):
	#	for j in self.gasses.keys:
	#		self.gasses[j].diffusion = self.catalyst_diffusion*(mp.sqrt(ref_mass*temperature)/mp.sqrt(ref_T*self.gasses[j].mass))
	#	for j in self.inertGasses.keys:
	#		pass
			
	#@property
	#def diffusion_matrix(self):
	#	return self._diffusion_matrix
	#@diffusion_matrix.setter
	#def diffusion_matrix(self, value):
	#	_diffusion_matrix = np.empty((len(self.gasses),3))
#
#		def diff_func(self,ref_mass,ref_T,mol_mass,ref_r):
#			return self.catalyst_diffusion*(mp.sqrt(ref_mass*temperature)/mp.sqrt(ref_T*mol_masss))
#
#		for k,j in enumerate(reac_input['Mass List'].split(',')):
#			pass
#
#	def reactor_diffusion_matrix(self, temperature):
#		D = np.empty((len(reac_input['Mass List'].split(',')),3))
#			
#		Dout = []
#		Din = []
#		constantTemp = reac_input['Reactor Temperature']
#		constantTemp = Constant(constantTemp)
#		testTemp = reac_input['Reactor Temperature']
#		testTemp = Constant(testTemp)
#
#
#		compMass_list = reac_input['Mass List'].split(',')
#	
#		for k,j in enumerate(reac_input['Mass List'].split(',')):
#			for k_2,j_2 in enumerate(ref_rate):
#				D[k,k_2] = diff_func(reac_input['Reference Mass'],reac_input['Reference Temperature'],float(j),j_2) ###??? should this (np.sum(r_param)**2) be here? This puts it in dimensional form!
#				if k_2 == 0:
#					Dout.append(Constant(D[k,k_2]))
#				if k_2 == 1:
#					Din.append(Constant(D[k,k_2]))

def initializeVariableDictionaries(nec,moles,V_nu,u_nu,un_nu):
	
	"""For all monitored parameters (gasses/surface species), establish necessary test and trial functions, as well as dictionaries for value storing"""

	graph_data = {}
	v_d = {}
	u_d = {}
	u_nd = {}
	surf_data = {}
	sens_data = {}	
	cat_data = {}

	species_count = nec['gas_num']+nec['surf_num']

	for k in range(0,nec['molecules_in_gas_phase']):
		graph_data['conVtime_'+str(k)] = []
		cat_data['conVtime_'+str(k)] = []
		sens_data['conVtime_'+str(k)] = []
	
	graph_data['conVtime_'+str(species_count)] = []
	sens_data['conVtime_'+str(species_count)] = []
	
	for kj in range(nec['molecules_in_gas_phase'],len(nec['reactants'])):
		surf_data['conVtime_'+str(kj)] = []
	
	tempA = TestFunctions(V_nu)
	tempB = split(u_nu)
	tempC = split(un_nu)
	
	for kit in range(0,int(moles)):
		v_d['v_'+str(kit+1)] = tempA[kit]
		u_d['u_'+str(kit+1)] = tempB[kit]
		u_nd['u_n'+str(kit+1)] = tempC[kit]
	
	return graph_data,v_d,u_d,u_nd,sens_data,surf_data,cat_data

#def fluxGeneration(reactor,gasses,reactants,pulse_size,Diff,voidage,dx,radius,dx2_r,outScale):
#	
#	"""Scale the output values of flux with the constant found here (messy now due to different trials)"""
#
#	to_flux = []
#
#	if reactor == 'tap':
#
#		for k in range(0,gasses):
#			if outScale.lower() == 'true':
#				to_flux.append( 2*(Diff[k][0] /(dx)) * (radius**2)*3.14159 / pulse_size ) #0.53*   ## 
#			else:
#				to_flux.append(2*(Diff[k][0] /(dx)) * (radius**2)*3.14159 )
#	elif reactor == 't_pfr' or 't_pfr_diff':
#		print('to flux generation')
#
#		for k in range(0,gasses):
#			to_flux.append(1)
#		
#		to_flux.append(1)
#	
#	return to_flux
#
#def pulse_functions(reactions_n,number_of_inerts):
#
#	rate_array, reactions, reactants, rev_irr = variational_list_parsing(reactions_n)
#	
#	gas_num = 0
#
#	gas_num = len(reactants)
#	molecules_in_gas_phase = 0
#
#	for k in reactants:
#		if '*' in k:
#			break
#		else:
#			molecules_in_gas_phase += 1 
#	
#	gas_molecules = reactants[:gas_num]
#	
#	gas_num = len(reactants)+int(number_of_inerts)
#
#	tail = len(reactants)
#
#
#	Fnew = '+ '	
#
#	for k in range(0,molecules_in_gas_phase):
#		Fnew += "-intensConst['inten_"+str(k)+"']*b0Test2*exp(-(constantT - timeConst['sST_"+str(k)+"'])*(constantT - timeConst['sST_"+str(k)+"'])/(4*0.00000000001))*v_d['v_"+str(k+1)+"']*dx"
#		if k < molecules_in_gas_phase-1:
#			Fnew += ' + '
#		elif int(number_of_inerts) > 0:
#			Fnew += ' + '
#	for k in range(0,int(number_of_inerts)):
#		Fnew += "-intensConst['inten_"+str(molecules_in_gas_phase+k)+"']*b0Test2*exp(-(constantT - timeConst['sST_"+str(molecules_in_gas_phase+k)+"'])*(constantT - timeConst['sST_"+str(molecules_in_gas_phase+k)+"'])/(4*0.00000000001))*v_d['v_"+str(len(reactants)+1+k)+"']*dx"
#		if k < int(number_of_inerts)-1:
#			Fnew += ' + '
#
#	return Fnew
#
#def batch_functions(reactions_n):
#
#	rate_array, reactions, reactants, rev_irr = variational_list_parsing(reactions_n)
#	
#	gas_num = 0
#
#	gas_num = len(reactants)
#	molecules_in_gas_phase = 0
#
#	for k in reactants:
#		if '*' in k:
#			break
#		else:
#			molecules_in_gas_phase += 1 
#	
#	gas_molecules = reactants[:gas_num]
#	
#	gas_num = len(reactants)
#
#	tail = len(reactants)
#
#
#	Fnew = '+ '	
#
#	for k in range(0,molecules_in_gas_phase):
#		Fnew += "-intensConst['inten_"+str(k)+"']*exp(-(constantT - timeConst['sST_"+str(k)+"'])*(constantT - timeConst['sST_"+str(k)+"'])/(4*0.00000000001))*v_d['v_"+str(k+1)+"']*dx"
#		if k < molecules_in_gas_phase-1:
#			Fnew += ' + '
#
#	return Fnew
#
#def surface_functions(reactions_n,number_of_inerts):
#	
#	rate_array, reactions, reactants, rev_irr = variational_list_parsing(reactions_n)
#	
#	gas_num = 0
#
#	gas_num = len(reactants)
#	molecules_in_gas_phase = 0
#
#	for k in reactants:
#		if '*' in k:
#			break
#		else:
#			molecules_in_gas_phase += 1 
#	
#	gas_molecules = reactants[:gas_num]
#	
#	gas_num = len(reactants)+int(number_of_inerts)
#
#	tail = len(reactants)
#
#
#	Fnew = ' '	
#
#	for k in range(molecules_in_gas_phase,len(reactants)):
#		Fnew += "- inComp["+str(k-molecules_in_gas_phase)+"]*exp(-(constantT-0.001)*(constantT-0.001)/(4*0.000000001))*v_d['v_"+str(1+k)+"']*dx "
#
#	return Fnew
#
#def batch_surface_functions(reactions_n):
#	
#	rate_array, reactions, reactants, rev_irr = variational_list_parsing(reactions_n)
#	
#	gas_num = 0
#
#	gas_num = len(reactants)
#	molecules_in_gas_phase = 0
#
#	for k in reactants:
#		if '*' in k:
#			break
#		else:
#			molecules_in_gas_phase += 1 
#	
#	gas_molecules = reactants[:gas_num]
#	
#	gas_num = len(reactants)
#
#	tail = len(reactants)
#
#
#	Fnew = ' '	
#
#	for k in range(molecules_in_gas_phase,len(reactants)):
#		Fnew += "- inComp["+str(k-molecules_in_gas_phase)+"]*exp(-(constantT-0.001)*(constantT-0.001)/(4*0.000000001))*v_d['v_"+str(1+k)+"']*dx "
#
#	return Fnew
#################################
#		# Read the initial composition of the catalyst
#		if ',' in str(reac_input['Initial Surface Composition']):
#			rangesurface_species = list(reversed(reac_input['Initial Surface Composition'].split(',')))
#			reac_input['Initial Surface Composition'] = list(map(float, rangesurface_species))
#		else:
#			rangesurface_species = reac_input['Initial Surface Composition']
#		
#
#		#############################################################
#		############# PARSE INLET PULSE COMPOSITION #################
#		#############################################################
#		
#		if '/' in str(reac_input['Pulse Size']):
#			list_of_feed = reac_input['Pulse Size'].split('/')
#			list_of_time = reac_input['Pulse Time'].split('/')
#			pulse_variation = len(list_of_feed)
#			pulse_time = len(list_of_time)
#			reactant_feed = []
#			reactant_time = []
#			for k in range(0,len(reac_input['Pulse Size'].split('/'))):
#				reactant_feed.append(list_of_feed[k].split(','))
#				reactant_time.append(list_of_time[k].split(','))
#		else:
#			pulse_variation = 1
#			reactant_feed = []
#			reactant_time = []
#			
#			try:
#				#reactant_feed.append(reac_input['Pulse Size'].split(','))
#				reactant_feed = list(map(float, reac_input['Pulse Size'].split(',')))
#				reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
#			except AttributeError:
#				#reactant_feed.append(reac_input['Pulse Size'])
#				reactant_feed = list(map(float, reac_input['Pulse Size'].split(',')))
#				reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
