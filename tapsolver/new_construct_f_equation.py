
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

from fenics import *
import pandas as pd
import numpy as np

import sys
#from structures import TAPobject
from TAPobject import TAPobject


def new_construct_f_equation(TAPobject_data: TAPobject):

	"""
	# TAPobject_data.species.temperature["+str(temperature_number)+"])
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
		return "((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']))*v_d['v_"+species_name+"']*dx      + dk*TAPobject_data.diffusion_inclusion*("+"(zone_diffusions*sqrt(TAPobject_data.species.reference_mass*TAPobject_data.species.temperature)/sqrt(TAPobject_data.species.reference_temperature*TAPobject_data.species.gasses['"+species_name+"'].mass))"+"/(reference_voids*(TAPobject_data.reactor.length**2)) ) *dot(grad(theta*u_d['u_"+species_name+"']+(1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx"
	def add_advection(species_number):
		"""Add advective term for reactive or inert gas species introduced in the reactor"""
		return " (dk)*Constant(theta)*Constant(TAPobject_data.reactor.length)*div(advMulti*advTerm*u_d['u_"+species_name+"'])*v_d['v_"+species_name+"']*dx + " + " (dk)*Constant(1-theta)*Constant(np.sum(TAPobject_data.reactor.length))*div(advMulti*advTerm*u_nd['u_n"+species_name+"'])*v_d['v_"+species_name+"']*dx "
	# For each reactive gas species, add the diffusion and advection terms (when appropriate)
	for knum, k in enumerate(list(TAPobject_data.species.gasses.keys())):
		if knum < len(TAPobject_data.species.gasses)-1:
			F = F+add_diffusion(k)+' + '
			if TAPobject_data.species.advection != 0:
				F = F+add_advection(k)+' + '	
		else:
			F = F+add_diffusion(k)			
			if TAPobject_data.species.advection != 0:
				F = F+' + '+add_advection(k)

	# For each adspecies, define the time variation term 
	for k in list(TAPobject_data.species.adspecies.keys()):
		F = F + " + ((u_d['u_"+k+"'] - u_nd['u_n"+k+"']))*v_d['v_"+k+"']*dx "
		
	def make_g(process,direction,scale_magnitude):
		"""Add free energy reaction term for the elementary process specified"""
		if direction == 'f': # standard_parameters['Av']*
			return "((standard_parameters['kbt']*TAPobject_data.species.temperature*(1/TAPobject_data.parameter_scale**("+str(scale_magnitude)+")))/(standard_parameters['h']))*exp(-TAPobject_data.mechanism.processes["+str(process)+"].f.Ga/(standard_parameters['kbt']*TAPobject_data.species.temperature))"
		else:
			return "((standard_parameters['kbt']*TAPobject_data.species.temperature*(1/TAPobject_data.parameter_scale**("+str(scale_magnitude)+")))/standard_parameters['h'])*exp(-(TAPobject_data.mechanism.processes["+str(process)+"].f.Ga - TAPobject_data.mechanism.processes["+str(process)+"].f.dG)/(standard_parameters['kbt']*TAPobject_data.species.temperature))"
			
	def make_arr(process,direction):
		"""Add activation/Arrhenius based reaction term for the elementary process specified"""
		if direction == 'f':
			return "TAPobject_data.mechanism.processes["+str(process)+"].f.Ao*exp(-TAPobject_data.mechanism.processes["+str(process)+"].f.Ea/(standard_parameters['Rgas']*TAPobject_data.species.temperature))"
		else:
			return "TAPobject_data.mechanism.processes["+str(process)+"].b.Ao*exp(-TAPobject_data.mechanism.processes["+str(process)+"].b.Ea/(standard_parameters['Rgas']*TAPobject_data.species.temperature))"

	def make_constant(process,direction):
		"""Add rate constant for the elementary process specified"""
		if direction == 'f':
			return "TAPobject_data.mechanism.processes["+str(process)+"].f.k"
		else:
			return "TAPobject_data.mechanism.processes["+str(process)+"].b.k"

	def make_link(process):
		"""Add rate constant for the elementary process specified"""
		return "TAPobject_data.mechanism.kinetic_links["+str(process)+"]"

	# Read through the stoichiometric matrix (i.e. rate_array) and define the associated system of odes for the mechanism
	if TAPobject_data.mechanism.reactants != []:
		for k,z in enumerate(TAPobject_data.mechanism.matrix):
			irr = None	
			neg = []
			val_neg = []
			pos = []
			val_pos = []
			for j,v in enumerate(z):
				if v < 0:
					neg.append(TAPobject_data.mechanism.reactants[j])
					val_neg.append(v)
				elif v > 0:
					pos.append(TAPobject_data.mechanism.reactants[j])
					val_pos.append(v)

			together = neg+pos
		
			if TAPobject_data.mechanism.processes[k].f.use == 'G':
				new_neg = make_g(k,'f',abs(sum(val_neg))-1)
			elif TAPobject_data.mechanism.processes[k].f.use == 'E':
				new_neg = make_arr(k,'f')
			elif TAPobject_data.mechanism.processes[k].f.use == 'k':
				new_neg = make_constant(k,'f')
			elif TAPobject_data.mechanism.processes[k].f.use == 'link':
				new_neg = make_link(TAPobject_data.mechanism.processes[k].f.link)
			
			#sys.exit()
			for j,v in enumerate(neg): # "(parameter_scale**(-1.0))*"+ "(1/120)*"+
				new_neg = new_neg+"*(u_d['u_"+v+"']**"+str(abs(val_neg[j]))+")"
			
			if TAPobject_data.mechanism.processes[k].f.use == 'G':
				new_pos = make_g(k,'b',abs(sum(val_pos))-1)
			elif TAPobject_data.mechanism.processes[k].b.use == 'E':
				new_pos = make_arr(k,'b')
			elif TAPobject_data.mechanism.processes[k].b.use == 'k':
				new_pos = make_constant(k,'b')
			elif TAPobject_data.mechanism.processes[k].b.use == 'link':
				new_pos = make_link(TAPobject_data.mechanism.processes[k].b.link)
			else:
				irr = True

			for j,v in enumerate(pos):
				new_pos = new_pos+"*(u_d['u_"+v+"']**"+str(abs(val_pos[j]))+")"#
			
			for j,v in enumerate(together):
				if j < len(neg):
					if irr == None:
						F = F+"- dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+v+"']*dx)"+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+v+"']*dx)"
					else:
						F = F+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"
				else:
					if irr == None:
						F = F+"+ dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+v+"']*dx)"+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+v+"']*dx)"
					else:
						F = F+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+v+"']*dx)"

	if TAPobject_data.mechanism.reactants != []: # 
		for k,z in enumerate(TAPobject_data.mechanism.matrix):
			neg = []
			val_neg = []
			pos = []
			val_pos = []
			for j,v in enumerate(z):
				if v < 0:
					neg.append(TAPobject_data.mechanism.reactants[j])
					val_neg.append(v)
				elif v > 0:
					pos.append(TAPobject_data.mechanism.reactants[j])
					val_pos.append(v)
			
			together = neg+pos

			if TAPobject_data.mechanism.processes[k].f.use == 'G':
				new_neg = make_g(k,'f',abs(sum(val_neg))-1)
			elif TAPobject_data.mechanism.processes[k].f.use == 'E':
				new_neg = make_arr(k,'f')
			elif TAPobject_data.mechanism.processes[k].f.use == 'k':
				new_neg = make_constant(k,'f')
			elif TAPobject_data.mechanism.processes[k].f.use == 'link':
				new_neg = make_link(TAPobject_data.mechanism.processes[k].f.link)

			for j,v in enumerate(neg):
				new_neg = new_neg+"*(u_nd['u_n"+v+"']**"+str(abs(val_neg[j]))+")"
			
			if TAPobject_data.mechanism.processes[k].f.use == 'G':
				new_pos = make_g(k,'b',abs(sum(val_pos))-1)
			elif TAPobject_data.mechanism.processes[k].f.use == 'E':
				new_pos = make_arr(k,'b')
			elif TAPobject_data.mechanism.processes[k].f.use == 'k':
				new_pos = make_constant(k,'b')
			elif TAPobject_data.mechanism.processes[k].f.use == 'link':
				new_pos = make_link(TAPobject_data.mechanism.processes[k].f.link)

			else:
				irr = True
					
			for j,v in enumerate(pos):
				new_pos = new_pos+"*(u_nd['u_n"+v+"']**"+str(abs(val_pos[j]))+")"

			for j,v in enumerate(together):
				if j < len(neg):
					if irr == None:
						F = F+"- dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+v+"']*dx)"+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+v+"']*dx)"
					else:
						F = F+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+v+"']*dx)"
				else:
					if irr == None:
						F = F+"+ dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+v+"']*dx)"+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+v+"']*dx)"
					else:
						F = F+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+v+"']*dx)"
		
	def add_inert_diffusion(species_name):
		"""Add diffusion term for inert gas species introduced in the reactor"""		
		return "((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']))*v_d['v_"+species_name+"']*dx      + dk*("+"(zone_diffusions*sqrt(TAPobject_data.species.reference_mass*TAPobject_data.species.temperature)/sqrt(TAPobject_data.species.reference_temperature*TAPobject_data.species.inert_gasses['"+species_name+"'].mass))"+"/(reference_voids*TAPobject_data.reactor.length**2) ) *dot(grad(theta*u_d['u_"+species_name+"']+(1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx    " 
	
	# For each inert gas species, add the diffusion and advection terms (when appropriate)
	for knum, k in enumerate(list(TAPobject_data.species.inert_gasses.keys())):
		if knum < len(TAPobject_data.species.inert_gasses)-1:
			F = F+" +  "+add_inert_diffusion(k)+' + '
			if TAPobject_data.species.advection != 0:
				F = F+" +  "+add_advection(k)+' + '	
		else:
			F = F+" +  "+add_inert_diffusion(k)			
			if TAPobject_data.species.advection != 0:
				F = F+' + '+add_advection(k)
	#sys.exit()
	return F