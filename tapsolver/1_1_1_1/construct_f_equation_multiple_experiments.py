
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

from fenics import *
import pandas as pd
import numpy as np

import sys
#from structures import TAPobject
from TAPobject import TAPobject 
from mechanism_reactants import mechanism_reactants


def construct_f_equation_multiple_experiments(TAPobject_data: TAPobject):

	"""

	This function constructs the variational form for TAP simulations in FEniCS, 
	incorporating the time steps, Knudsen diffusion, advection, and elementary
	reactions for the reactive gasses, inert gasses, and catalyst adspecies.

	Args:

		TAPobject_data (TAPobject): An object of the TAP class that includes the
		reactor, mechanism, and reactor species details.
	
	Returns:

		F (str): A string defining the variational form of the simulation, which
		is subsequently evaluated in the forward or inverse problem scripts. 

	Implementor:

		Adam Yonge

	Link:
		
		https://en.wikipedia.org/wiki/Weak_formulation
	
	"""

	F = ''
	def add_diffusion(species_name):
		"""Add diffusion term for reactive gas species introduced in the reactor"""
		return "((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']))*v_d['v_"+species_name+"']*dx(0)      + dk*("+"(TAPobject_data.reactor_species.inert_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature["+str(TAPobject_data.reactor_species.gasses[species_name].temperature_used)+"])/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.gasses['"+species_name+"'].mass))"+"/(TAPobject_data.reactor.zone_voids[0]*(TAPobject_data.reactor.total_length**2)) ) *dot(grad(theta*u_d['u_"+species_name+"']+(1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx(0)     + ((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']) )*v_d['v_"+species_name+"']*dx(1)       + dk*("+"(TAPobject_data.reactor_species.catalyst_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature["+str(TAPobject_data.reactor_species.gasses[species_name].temperature_used)+"])/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.gasses['"+species_name+"'].mass))"+"/(TAPobject_data.reactor.zone_voids[1]*(TAPobject_data.reactor.total_length**2)) )*dot(grad(theta*u_d['u_"+species_name+"'] + (1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx(1) " 
	
	def add_advection(species_number):
		"""Add advective term for reactive or inert gas species introduced in the reactor"""
		return " (dk)*Constant(theta)*Constant(TAPobject_data.reactor.total_length)*div(advMulti*advTerm*u_d['u_"+species_name+"'])*v_d['v_"+species_name+"']*dx(0) + " + " (dk)*Constant(1-theta)*Constant(np.sum(TAPobject_data.reactor.total_length))*div(advMulti*advTerm*u_nd['u_n"+species_name+"'])*v_d['v_"+species_name+"']*dx(0) + " + " (dk)*Constant(theta)*Constant(TAPobject_data.reactor.total_length)*div(advMulti*advTerm*u_d['u_"+species_name+"'])*v_d['v_"+species_name+"']*dx(1) + " + " (dk)*Constant(1-theta)*Constant(TAPobject_data.reactor.total_length)*div(advMulti*advTerm*u_nd['u_n"+species_name+"'])*v_d['v_"+species_name+"']*dx(1)"
	
	# For each reactive gas species, add the diffusion and advection terms (when appropriate)
	for knum, k in enumerate(list(TAPobject_data.reactor_species.gasses.keys())):
		if knum < len(TAPobject_data.reactor_species.gasses)-1:
			F = F+add_diffusion(k)+' + '
			if TAPobject_data.reactor_species.advection != 0:
				F = F+add_advection(k)+' + '	
		else:
			F = F+add_diffusion(k)			
			if TAPobject_data.reactor_species.advection != 0:
				F = F+' + '+add_advection(k)
	##################3	
	# For each adspecies, define the time variation term 
	for k in list(TAPobject_data.reactor_species.adspecies.keys()):
		F = F + " + ((u_d['u_"+k+"'] - u_nd['u_n"+k+"']))*v_d['v_"+k+"']*dx(0) + ((u_d['u_"+k+"'] - u_nd['u_n"+k+"']) )*v_d['v_"+k+"']*dx(1)"

	def make_g(elementary_process,direction,scale_magnitude):
		"""Add free energy reaction term for the elementary process specified"""
		if direction == 'f':
			return "((standard_parameters['kbt']*TAPobject_data.reactor_species.temperature["+str(temperature_number)+"]*(1/TAPobject_data.parameter_scale**("+str(scale_magnitude)+")))/(standard_parameters['h']))*exp(-TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.Ga/(standard_parameters['kbt']*TAPobject_data.reactor_species.temperature["+str(TAPobject_data.reactor_species.gasses[species_name].temperature_used)+"]))"
		else:
			return "((standard_parameters['kbt']*TAPobject_data.reactor_species.temperature["+str(temperature_number)+"]*(1/TAPobject_data.parameter_scale**("+str(scale_magnitude)+")))/(standard_parameters['h']))*exp(-(TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.Ga - TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.dG)/(standard_parameters['kbt']*TAPobject_data.reactor_species.temperature["+str(TAPobject_data.reactor_species.gasses[species_name].temperature_used)+"]))"
			
	def make_arr(elementary_process,direction):
		"""Add activation/Arrhenius based reaction term for the elementary process specified"""
		if direction == 'f':
			return "TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.Ao*exp(-TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.Ea/(standard_parameters['Rgas']*TAPobject_data.reactor_species.temperature["+str(TAPobject_data.reactor_species.gasses[species_name].temperature_used)+"]))"
		else:
			return "TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].backward.Ao*exp(-TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].backward.Ea/(standard_parameters['Rgas']*TAPobject_data.reactor_species.temperature["+str(TAPobject_data.reactor_species.gasses[species_name].temperature_used)+"]))"

	def make_constant(elementary_process,direction):
		"""Add rate constant for the elementary process specified"""
		if direction == 'f':
			return "TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.k"
		else:
			return "TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].backward.k"

	def make_link(elementary_process):
		"""Add rate constant for the elementary process specified"""
		return "TAPobject_data.mechanism.kinetic_links["+str(elementary_process)+"]"
		
	# Read through the stoichiometric matrix (i.e. rate_array) and define the associated system of odes for the mechanism
	
	#if TAPobject_data.mechanism.reactants[temperature_number] != []:
	for temperature_number in TAPobject_data.reactor_species.temperature:
		for k,z in enumerate(TAPobject_data.mechanism.rate_array):
			irr = None	
			neg = []
			val_neg = []
			pos = []
			val_pos = []
			for j,v in enumerate(z):
				if v < 0:
					neg.append(TAPobject_data.mechanism.reactants[temperature_number][j])
					val_neg.append(v)
				elif v > 0:
					pos.append(TAPobject_data.mechanism.reactants[temperature_number][j])
					val_pos.append(v)

			together = neg+pos
			
			if TAPobject_data.mechanism.elementary_processes[k].forward.use == 'G':
				new_neg = make_g(k,'f',abs(sum(val_neg))-1)
			elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'E':
				new_neg = make_arr(k,'f')
			elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'k':
				new_neg = make_constant(k,'f')
			elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'link':
				new_neg = make_link(TAPobject_data.mechanism.elementary_processes[k].forward.link)
				
			for j,v in enumerate(neg):
				new_neg = new_neg+"*(u_d['u_"+v+"']**"+str(abs(val_neg[j]))+")"

			if TAPobject_data.mechanism.elementary_processes[k].backward.use == 'G':
				new_pos = make_g(k,'b',abs(sum(val_pos))-1)
			elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'E':
				new_pos = make_arr(k,'b')
			elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'k':
				new_pos = make_constant(k,'b')
			elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'link':
				new_pos = make_link(TAPobject_data.mechanism.elementary_processes[k].backward.link)
			else:
				irr = True
			for j,v in enumerate(pos):
				new_pos = new_pos+"*(u_d['u_"+v+"']**"+str(abs(val_pos[j]))+")"#

			for j,v in enumerate(together):
				if j < len(neg):
					if irr == None:
						F = F+"- dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+v+"']*dx(1))"+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+v+"']*dx(1))"
					else:
						F = F+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
				else:
					if irr == None:
						F = F+"+ dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+v+"']*dx(1))"+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+v+"']*dx(1))"
					else:
						F = F+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+v+"']*dx(1))"

		#if TAPobject_data.mechanism.reactants[temperature_number] != []:
			for k,z in enumerate(TAPobject_data.mechanism.rate_array):
				neg = []
				val_neg = []
				pos = []
				val_pos = []
				for j,v in enumerate(z):
					if v < 0:
						neg.append(TAPobject_data.mechanism.reactants[temperature_number][j])
						val_neg.append(v)
					elif v > 0:
						pos.append(TAPobject_data.mechanism.reactants[temperature_number][j])
						val_pos.append(v)
				
				together = neg+pos

				if TAPobject_data.mechanism.elementary_processes[k].forward.use == 'G':
					new_neg = make_g(k,'f',abs(sum(val_neg))-1)
				elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'E':
					new_neg = make_arr(k,'f')
				#elif TAPobject_data.mechanism.elementary_processes[k].forward.link['variable'] != None:
				#	new_neg = make_link(k,'f')
				elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'k':
					new_neg = make_constant(k,'f')
				elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'link':
					new_neg = make_link(TAPobject_data.mechanism.elementary_processes[k].forward.link)

				for j,v in enumerate(neg):
					new_neg = new_neg+"*(u_nd['u_n"+v+"']**"+str(abs(val_neg[j]))+")"
				
				if TAPobject_data.mechanism.elementary_processes[k].backward.use == 'G':
					new_pos = make_g(k,'b',abs(sum(val_pos))-1)
				elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'E':
					new_pos = make_arr(k,'b')
				#elif TAPobject_data.mechanism.elementary_processes[k].backward.link['variable'] != 0:
				#	new_pos = make_link(k,'b')
				elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'k':
					new_pos = make_constant(k,'b')
				elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'link':
					new_pos = make_link(TAPobject_data.mechanism.elementary_processes[k].backward.link)
				else:
					irr = True
					
				for j,v in enumerate(pos):
					new_pos = new_pos+"*(u_nd['u_n"+v+"']**"+str(abs(val_pos[j]))+")"

				for j,v in enumerate(together):
					if j < len(neg):
						if irr == None:
							F = F+"- dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+v+"']*dx(1))"+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+v+"']*dx(1))"
						else:
							F = F+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+v+"']*dx(1))"
					else:
						if irr == None:
							F = F+"+ dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+v+"']*dx(1))"+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+v+"']*dx(1))"
						else:
							F = F+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+v+"']*dx(1))"
		
		def add_inert_diffusion(species_name):
			"""Add diffusion term for inert gas species introduced in the reactor"""
			return "((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']))*v_d['v_"+species_name+"']*dx(0)      + dk*("+"(TAPobject_data.reactor_species.inert_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature["+str(TAPobject_data.reactor_species.inert_gasses[species_name].temperature_used)+"])/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.inert_gasses['"+species_name+"'].mass))"+"/(TAPobject_data.reactor.zone_voids[0]*TAPobject_data.reactor.total_length**2) ) *dot(grad(theta*u_d['u_"+species_name+"']+(1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx(0)     + ((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']) )*v_d['v_"+species_name+"']*dx(1)       + dk*("+"(TAPobject_data.reactor_species.catalyst_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature["+str(TAPobject_data.reactor_species.inert_gasses[species_name].temperature_used)+"])/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.inert_gasses['"+species_name+"'].mass))"+"/(TAPobject_data.reactor.zone_voids[1]*TAPobject_data.reactor.total_length**2) )*dot(grad(theta*u_d['u_"+species_name+"'] + (1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx(1) " 
		
		# For each inert gas species, add the diffusion and advection terms (when appropriate)
		for knum, k in enumerate(list(TAPobject_data.reactor_species.inert_gasses.keys())):
			if knum < len(TAPobject_data.reactor_species.inert_gasses)-1:
				F = F+" +  "+add_inert_diffusion(k)+' + '
				if TAPobject_data.reactor_species.advection != 0:
					F = F+" +  "+add_advection(k)+' + '	
			else:
				F = F+" +  "+add_inert_diffusion(k)			
				if TAPobject_data.reactor_species.advection != 0:
					F = F+' + '+add_advection(k)
	return F