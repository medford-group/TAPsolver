
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
import sys
from structures import TAPobject
from .mechanism_reactants import mechanism_reactants


def construct_f_equation(TAPobject_data: TAPobject):

	active_site_names = ['*','^','@','#']

	F = ''
	rate_array = TAPobject_data.mechanism.rate_array
	reactions = TAPobject_data.mechanism.reactions
	reactants = TAPobject_data.mechanism.reactants
	number_of_inerts = len(TAPobject_data.reactor_species.inert_gasses)

	gas_num = 0

	gas_num = len(reactants)
	molecules_in_gas_phase = 0

	for k in reactants:
		if '*' in k:
			break
		else:
			molecules_in_gas_phase += 1 
	
	gas_molecules = reactants[:gas_num]
	
	gas_num = len(reactants)+int(number_of_inerts)

	tail = len(reactants)

	# TAPobject_data.reactor.zone_voids[0]
	# reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))
	# "+"(TAPobject_data.reactor_species.inert_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature)/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.gasses["+species_name+"].mass))"+"
	# "+"(TAPobject_data.reactor_species.catalyst_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature)/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.gasses["+species_name+"].mass))"+"


	def add_diffusion(species_name):
		return "((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']))*v_d['v_"+species_name+"']*dx(0)      + dk*("+"(TAPobject_data.reactor_species.inert_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature)/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.gasses['"+species_name+"'].mass))"+"/(TAPobject_data.reactor.zone_voids[0]*TAPobject_data.reactor.total_length**2) ) *dot(grad(theta*u_d['u_"+species_name+"']+(1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx(0)     + ((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']) )*v_d['v_"+species_name+"']*dx(1)       + dk*("+"(TAPobject_data.reactor_species.catalyst_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature)/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.gasses['"+species_name+"'].mass))"+"/(TAPobject_data.reactor.zone_voids[1]*TAPobject_data.reactor.total_length**2) )*dot(grad(theta*u_d['u_"+species_name+"'] + (1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx(1) " 
	
	def add_advection(species_number):
		return " (dk)*Constant(theta)*Constant(TAPobject_data.reactor.total_length)*div(advMulti*advTerm*u_d['u_"+species_name+"'])*v_d['v_"+species_name+"']*dx(0) + " + " (dk)*Constant(1-theta)*Constant(np.sum(TAPobject_data.reactor.total_length))*div(advMulti*advTerm*u_nd['u_n"+species_name+"'])*v_d['v_"+species_name+"']*dx(0) + " + " (dk)*Constant(theta)*Constant(TAPobject_data.reactor.total_length)*div(advMulti*advTerm*u_d['u_"+species_name+"'])*v_d['v_"+species_name+"']*dx(1) + " + " (dk)*Constant(1-theta)*Constant(TAPobject_data.reactor.total_length)*div(advMulti*advTerm*u_nd['u_n"+species_name+"'])*v_d['v_"+species_name+"']*dx(1)"
	
	for knum, k in enumerate(list(TAPobject_data.reactor_species.gasses.keys())):
		if knum < len(TAPobject_data.reactor_species.gasses)-1:
			F = F+add_diffusion(k)+' + '
			if TAPobject_data.reactor_species.advection != 0:
				F = F+add_advection(k)+' + '	
		else:
			F = F+add_diffusion(k)			
			if TAPobject_data.reactor_species.advection != 0:
				F = F+' + '+add_advection(k)

	for k in list(TAPobject_data.reactor_species.gasses.keys()):
		F = F + " + ((u_d['u_"+k+"'] - u_nd['u_n"+k+"']))*v_d['v_"+k+"']*dx(0) + ((u_d['u_"+k+"'] - u_nd['u_n"+k+"']) )*v_d['v_"+k+"']*dx(1)"
		
	

	def make_g(elementary_process,direction):
		if direction == 'f':
			return "standard_parameters['kbt']*TAPobject_data.reactor_species.temperature/standard_parameters['h']*exp(-TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.Ga['value']/(standard_parameters['kbt']*TAPobject_data.reactor_species.temperature))"
		else:
			return "standard_parameters['kbt']*TAPobject_data.reactor_species.temperature/standard_parameters['h']*exp(-(TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.Ga['value'] - TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.dG['value'])/(standard_parameters['kbt']*TAPobject_data.reactor_species.temperature))"
			
	def make_arr(elementary_process,direction):
		#standard_parameters['Rgas']
		if direction == 'f':
			return "TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.Ao['value']*exp(-TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.Ea['value']/(standard_parameters['Rgas']*TAPobject_data.reactor_species.temperature))"
		else:
			return "TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].backward.Ao['value']*exp(-TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].backward.Ea['value']/(standard_parameters['Rgas']*TAPobject_data.reactor_species.temperature))"
		
		#return 'r_Ao["Ao'+direction+str(parameter_number)+'"]*exp(-r_Ea["Ea'+direction+str(parameter_number)+'"]/(Rgas*constantTemp))'

	def make_link(elementary_process,direction):
		return 'r_links[kineticLinks["k'+direction+str(parameter_number)+'"]]'

	def make_constant(elementary_process,direction):
		if direction == 'f':
			return "TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].forward.k['value']"
		else:
			return "TAPobject_data.mechanism.elementary_processes["+str(elementary_process)+"].backward.k['value']"


	for k,z in enumerate(rate_array):
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
		
		if TAPobject_data.mechanism.elementary_processes[k].forward.use == 'G':
			#if k in gForward:#,r_Ga_in,r_dG_in
			new_neg = make_g(k,'f')
		elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'E':
			#elif k in arrForward:
			new_neg = make_arr(k,'f')
		#elif TAPobject_data.mechanism.elementary_processes[k].forward.link['variable'] != None:
		#	#elif k in linkForward:
		#	new_neg = make_link(k,'f')
		elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'k':
			new_neg = make_constant(k,'f')
				
		for j,v in enumerate(neg):
			new_neg = new_neg+"*(u_d['u_"+v+"']**"+str(abs(val_neg[j]))+")"

		if TAPobject_data.mechanism.elementary_processes[k].backward.use == 'G':
			new_pos = make_g(k,'b')
		elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'E':
			new_pos = make_arr(k,'b')
		#elif TAPobject_data.mechanism.elementary_processes[k].backward.link['variable'] != 0:
		#	new_pos = make_link(k,'b')
		elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'k':
			new_pos = make_constant(k,'b')
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


#######################################################

	for k,z in enumerate(rate_array):
				
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

		if TAPobject_data.mechanism.elementary_processes[k].forward.use == 'G':
			new_neg = make_g(k,'f')
		elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'E':
			new_neg = make_arr(k,'f')
		#elif TAPobject_data.mechanism.elementary_processes[k].forward.link['variable'] != None:
		#	new_neg = make_link(k,'f')
		elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'k':
			new_neg = make_constant(k,'f')

		for j,v in enumerate(neg):
			new_neg = new_neg+"*(u_nd['u_n"+v+"']**"+str(abs(val_neg[j]))+")"
			
		if TAPobject_data.mechanism.elementary_processes[k].backward.use == 'G':
			new_pos = make_g(k,'b')
		elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'E':
			new_pos = make_arr(k,'b')
		#elif TAPobject_data.mechanism.elementary_processes[k].backward.link['variable'] != 0:
		#	new_pos = make_link(k,'b')
		elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'k':
			new_pos = make_constant(k,'b')
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
	
#######################################################
	
		
	def add_inert_diffusion(species_name):
		return "((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']))*v_d['v_"+species_name+"']*dx(0)      + dk*("+"(TAPobject_data.reactor_species.inert_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature)/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.inert_gasses['"+species_name+"'].mass))"+"/(TAPobject_data.reactor.zone_voids[0]*TAPobject_data.reactor.total_length**2) ) *dot(grad(theta*u_d['u_"+species_name+"']+(1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx(0)     + ((u_d['u_"+species_name+"'] - u_nd['u_n"+species_name+"']) )*v_d['v_"+species_name+"']*dx(1)       + dk*("+"(TAPobject_data.reactor_species.catalyst_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature)/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.inert_gasses['"+species_name+"'].mass))"+"/(TAPobject_data.reactor.zone_voids[1]*TAPobject_data.reactor.total_length**2) )*dot(grad(theta*u_d['u_"+species_name+"'] + (1-theta)*u_nd['u_n"+species_name+"']), grad(v_d['v_"+species_name+"']))*dx(1) " 
	
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

	#for k in range(0,int(number_of_inerts)):
	#		
	#	F = F + "+  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']))*v_d['v_"+str(len(reactants)+1+k)+"']*dx(0)        +  Dout["+str(molecules_in_gas_phase+k)+"]* dk *( 1/Constant(eb[0]*np.sum(r_param)**2) ) * dot(grad( theta*u_d['u_"+str(len(reactants)+1+k)+"'] + (1-theta)*u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(0)     +  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']) )*v_d['v_"+str(len(reactants)+1+k)+"']*dx(1)             + dk * (Din["+str(molecules_in_gas_phase+k)+"]  /Constant(eb[1]*np.sum(r_param)**2) ) * dot(grad( theta*u_d['u_"+str(len(reactants)+1+k)+"'] + (1-theta)*u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(1)"
	#		
	#	if TAPobject_data.reactor_species.advection != 0:
	#		F = F+ "+ dk *Constant(theta)* Constant((np.sum(r_param))) * dot(advMulti*advTerm, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)        + dk *Constant(theta)* Constant((np.sum(r_param))) * dot(advMulti*advTerm, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(1) "

	
