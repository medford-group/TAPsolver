
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved



class elementary_process_details():
	
	# default (str): Define what kind of parameters should be used for the simulation or inverse problem ('k': rate constants, 'g': free energy, 'a': activation, 'l': linked)
	
	def __init__(self):
		self.k = None
		self.Ao = None
		self.Ea = None
		self.Ga = None
		self.dG = None
		self.link = None
		#self.link = {'variable':0}
		self.use = None # 'None' or 'k' # or 'E' or 'G' or 'l' or []
		self.lower_bound = None
		self.upper_bound = None