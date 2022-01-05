
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

from fenics import *
from fenics_adjoint import *
from .mechanism import mechanism
from .mechanism_construction import *
import copy

def load_kinetics(mechanism_data: mechanism):

	mechanism_copy = copy.deepcopy(mechanism_data)

	if mechanism_copy.rate_array == None:
		mechanism_copy = mechanism_constructor(mechanism_copy)

	mechanism_processes = mechanism_copy.elementary_processes.keys()
	
	r_Ga_in = {}
	r_dG_in = {}
	r_const = {}
	r_Ao = {}
	r_Ea = {}
	r_links = {}

	for j in mechanism_processes:
		if mechanism_copy.elementary_processes[j].forward.Ga['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'g':
			r_Ga_in["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.Ga['value'])
		if mechanism_copy.elementary_processes[j].forward.dG['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'g':
			r_dG_in["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.dG['value'])
		if mechanism_copy.elementary_processes[j].forward.k['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'k':
			r_const["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.k['value'])
		if mechanism_copy.elementary_processes[j].forward.Ao['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'a':
			r_Ao["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.Ao['value'])
		if mechanism_copy.elementary_processes[j].forward.Ea['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'a':
			r_Ea["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.Ea['value'])

		if mechanism_copy.elementary_processes[j].backward.k['value'] != None and mechanism_copy.elementary_processes[j].backward.default == 'k':
			r_const["kb"+str(j)] = Constant(mechanism_copy.elementary_processes[j].backward.k['value'])
		if mechanism_copy.elementary_processes[j].backward.Ao['value'] != None and mechanism_copy.elementary_processes[j].backward.default == 'a':
			r_Ao["kb"+str(j)] = Constant(mechanism_copy.elementary_processes[j].backward.Ao['value'])
		if mechanism_copy.elementary_processes[j].backward.Ea['value'] != None and mechanism_copy.elementary_processes[j].backward.default == 'a':
			r_Ea["kb"+str(j)] = Constant(mechanism_copy.elementary_processes[j].backward.Ea['value'])

	for jnum,j in enumerate(r_Ga_in):
		r_Ga_in[j] = Constant(r_Ga_in[j])
	
	for j in r_dG_in:
		r_dG_in[j] = Constant(r_dG_in[j])
	
	for j in r_const:
		r_const[j] = Constant(r_const[j])
	
	for j in r_Ao:
		r_Ao[j] = Constant(r_Ao[j])
	
	for j in r_Ea:
		r_Ea[j] = Constant(r_Ea[j])

	for j in r_links:
		r_links[j] = Constant(r_links[j])

	kinetic_information = {'r_Ga_in': r_Ga_in, 'r_dG_in': r_dG_in, 'r_const': r_const, 'r_Ao': r_Ao, 'r_Ea': r_Ea, 'r_links': r_links}

	return kinetic_information