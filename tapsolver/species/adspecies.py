
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

class adspecies():
	def __init__(self):
		self.conc = 0
		self.noise = 0 
		self.sigma = 1
		self.profile = 'uniform' # uniform # linear # bounded
		self.p_parms = {} # 1 # ax + b # [< + <]  
		self.zones = [1]