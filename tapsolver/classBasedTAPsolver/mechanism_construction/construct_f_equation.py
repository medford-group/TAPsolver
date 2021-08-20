
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
import sys
from structures import mechanism, initial_conditions, reactor
from .mechanism_reactants import mechanism_reactants


def construct_f_equation(mechanism_data: mechanism,initial_conditions_data: initial_conditions, reactor_data: reactor):

	active_site_names = ['*','^','@','#']

	F = ''
	rate_array = mechanism_data.rate_array
	reactions = mechanism_data.reactions
	reactants = mechanism_data.reactants
	number_of_inerts = len(initial_conditions_data.inertGasses.keys())

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

	def add_diffusion(species_number):
		return "((u_d['u_"+str(species_number+1)+"'] - u_nd['u_n"+str(species_number+1)+"']))*v_d['v_"+str(species_number+1)+"']*dx(0)      + dk*(Dout["+str(species_number)+"]/(eb[0]*np.sum(r_param)**2) ) *dot(grad(theta*u_d['u_"+str(species_number+1)+"']+(1-theta)*u_nd['u_n"+str(species_number+1)+"']), grad(v_d['v_"+str(species_number+1)+"']))*dx(0)     + ((u_d['u_"+str(species_number+1)+"'] - u_nd['u_n"+str(species_number+1)+"']) )*v_d['v_"+str(species_number+1)+"']*dx(1)       + dk*(Din["+str(species_number)+"]/(eb[1]*np.sum(r_param)**2) )*dot(grad(theta*u_d['u_"+str(species_number+1)+"'] + (1-theta)*u_nd['u_n"+str(species_number+1)+"']), grad(v_d['v_"+str(species_number+1)+"']))*dx(1) " 
	
	def add_advection(species_number):
		return " (dk)*Constant(theta)*Constant(np.sum(r_param))*div(advMulti*advTerm*u_d['u_"+str(species_number+1)+"'])*v_d['v_"+str(species_number+1)+"']*dx(0) + " + " (dk)*Constant(1-theta)*Constant(np.sum(r_param))*div(advMulti*advTerm*u_nd['u_n"+str(species_number+1)+"'])*v_d['v_"+str(species_number+1)+"']*dx(0) + " + " (dk)*Constant(theta)*Constant(np.sum(r_param))*div(advMulti*advTerm*u_d['u_"+str(species_number+1)+"'])*v_d['v_"+str(species_number+1)+"']*dx(1) + " + " (dk)*Constant(1-theta)*Constant(np.sum(r_param))*div(advMulti*advTerm*u_nd['u_n"+str(species_number+1)+"'])*v_d['v_"+str(species_number+1)+"']*dx(1)"

	
	for k in range(0,molecules_in_gas_phase):
		if k < (molecules_in_gas_phase-1):
			F = F+add_diffusion(k)+' + '
			if reactor_data.advection != 0:
				F = F+add_advection(k)+' + '	
		else:
			F = F+add_diffusion(k)			
			if reactor_data.advection != 0:
				F = F+' + '+add_advection(k)

	for k in range(molecules_in_gas_phase,len(reactants)):
		F = F + " + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)"
			
	def make_g(parameter_number,direction):
		if direction == 'f':
			return '(kbt*constantTemp/hb)*exp(-r_Ga_in["Ga'+str(parameter_number)+'"])'
		else:
			return '(kbt*constantTemp/hb)*exp((-(r_Ga_in["Ga'+str(parameter_number)+'"]-r_dG_in["dG'+str(parameter_number)+'"])))'

	def make_arr(parameter_number,direction):
		return 'r_Ao["Ao'+direction+str(parameter_number)+'"]*exp(-r_Ea["Ea'+direction+str(parameter_number)+'"]/(Rgas*constantTemp))'

	def make_link(parameter_number,direction):
		return 'r_links[kineticLinks["k'+direction+str(parameter_number)+'"]]'

	def make_constant(parameter_number,direction):
		return 'r_const["k'+direction+str(parameter_number)+'"]'


	for k,z in enumerate(rate_array):
		irr = None	
		neg = []
		val_neg = []
		pos = []
		val_pos = []
		for j,v in enumerate(z):
			if v < 0:
				neg.append(j)
				val_neg.append(v)
			elif v > 0:
				pos.append(j)
				val_pos.append(v)

		together = neg+pos
		
		if mechanism_data.elementaryProcesses[k].forward.Ga['value'] != None:
			#if k in gForward:#,r_Ga_in,r_dG_in
			new_neg = make_g(k,'f')
		elif mechanism_data.elementaryProcesses[k].forward.Ao['value'] != None:	
			#elif k in arrForward:
			new_neg = make_arr(k,'f')
		elif mechanism_data.elementaryProcesses[k].forward.link['variable'] != None:
			#elif k in linkForward:
			new_neg = make_link(k,'f')
		else:
			new_neg = make_constant(k,'f')
				
		for j,v in enumerate(neg):
			new_neg = new_neg+"*(u_d['u_"+str(v+1)+"']**"+str(abs(val_neg[j]))+")"

		if mechanism_data.elementaryProcesses[k].backward.Ga['value'] != None:
			new_pos = make_g(k,'b')
		elif mechanism_data.elementaryProcesses[k].backward.Ao['value'] != None:
			new_pos = make_arr(k,'b')
		elif mechanism_data.elementaryProcesses[k].backward.link['variable'] != 0:
			new_pos = make_link(k,'b')
		elif mechanism_data.elementaryProcesses[k].backward.k['value'] != None:
			new_pos = make_constant(k,'b')
		else:
			irr = True

		for j,v in enumerate(pos):
			new_pos = new_pos+"*(u_d['u_"+str(v+1)+"']**"+str(abs(val_pos[j]))+")"#

		for j,v in enumerate(together):
			if j < len(neg):
				if irr == None:
					F = F+"- dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
				else:
					F = F+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
			else:
				if irr == None:
					F = F+"+ dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
				else:
					F = F+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"

#######################################################

	for k,z in enumerate(rate_array):
				
		neg = []
		val_neg = []
		pos = []
		val_pos = []
		for j,v in enumerate(z):
			if v < 0:
				neg.append(j)
				val_neg.append(v)
			elif v > 0:
				pos.append(j)
				val_pos.append(v)
			
		together = neg+pos

		if mechanism_data.elementaryProcesses[k].forward.Ga['value'] != None:
			new_neg = make_g(k,'f')
		elif mechanism_data.elementaryProcesses[k].forward.Ao['value'] != None:	
			new_neg = make_arr(k,'f')
		elif mechanism_data.elementaryProcesses[k].forward.link['variable'] != None:	
			new_neg = make_link(k,'f')
		else:
			new_neg = make_constant(k,'f')

		for j,v in enumerate(neg):
			new_neg = new_neg+"*(u_nd['u_n"+str(v+1)+"']**"+str(abs(val_neg[j]))+")"
			
		if mechanism_data.elementaryProcesses[k].backward.Ga['value'] != None:
			new_pos = make_g(k,'b')
		elif mechanism_data.elementaryProcesses[k].backward.Ao['value'] != None:
			new_pos = make_arr(k,'b')
		elif mechanism_data.elementaryProcesses[k].backward.link['variable'] != 0:
			new_pos = make_link(k,'b')
		elif mechanism_data.elementaryProcesses[k].backward.k['value'] != None:
			new_pos = make_constant(k,'b')
		else:
			irr = True
				
		for j,v in enumerate(pos):
			new_pos = new_pos+"*(u_nd['u_n"+str(v+1)+"']**"+str(abs(val_pos[j]))+")"

		for j,v in enumerate(together):
			if j < len(neg):
				if irr == None:
					F = F+"- dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
				else:
					F = F+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
			else:
				if irr == None:
					F = F+"+ dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
				else:
					F = F+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
	
#######################################################

		for k in range(0,int(number_of_inerts)):
			
			F = F + "+  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']))*v_d['v_"+str(len(reactants)+1+k)+"']*dx(0)        +  Dout["+str(molecules_in_gas_phase+k)+"]* dk *( 1/Constant(eb[0]*np.sum(r_param)**2) ) * dot(grad( theta*u_d['u_"+str(len(reactants)+1+k)+"'] + (1-theta)*u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(0)     +  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']) )*v_d['v_"+str(len(reactants)+1+k)+"']*dx(1)             + dk * (Din["+str(molecules_in_gas_phase+k)+"]  /Constant(eb[1]*np.sum(r_param)**2) ) * dot(grad( theta*u_d['u_"+str(len(reactants)+1+k)+"'] + (1-theta)*u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(1)"
			
			if reactor_data.advection != 0:
				F = F+ "+ dk *Constant(theta)* Constant((np.sum(r_param))) * dot(advMulti*advTerm, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)        + dk *Constant(theta)* Constant((np.sum(r_param))) * dot(advMulti*advTerm, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(1) "

	
###################
	element = '['
	for k in range(0,len(reactants)+int(number_of_inerts)):
		if (k+1) < len(reactants)+int(number_of_inerts):
			element = element + 'P1,'
		else:
			element = element + 'P1]'

	return F, element