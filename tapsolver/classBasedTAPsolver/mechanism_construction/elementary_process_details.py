
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved



class elementary_process_details():
	
	# default (str): Define what kind of parameters should be used for the simulation or inverse problem ('k': rate constants, 'g': free energy, 'a': activation, 'l': linked)
	
	def __init__(self):
		self.k = {'value':None,'fd':'fixed'}
		self.Ao = {'value':None,'fd':'fixed'}
		self.Ea = {'value':None,'fd':'fixed'}
		self.Ga = {'value':None,'fd':'fixed'}
		self.dG = {'value':None,'fd':'fixed'}
		self.link = {'variable':0}
		self.default = 'k'