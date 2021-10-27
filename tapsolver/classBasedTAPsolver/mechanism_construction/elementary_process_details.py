
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved



class elementary_process_details():
	
	# default (str): Define what kind of parameters should be used for the simulation or inverse problem ('k': rate constants, 'g': free energy, 'a': activation, 'l': linked)
	
	def __init__(self):
		self.k = {'value':None,'fd':'fit'}
		self.Ao = {'value':None,'fd':'fix'}
		self.Ea = {'value':None,'fd':'fix'}
		self.Ga = {'value':None,'fd':'fix'}
		self.dG = {'value':None,'fd':'fix'}
		#self.link = {'variable':0}
		self.use = 'k' # or 'E' or 'G' or 'l' or []