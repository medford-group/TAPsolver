
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved


class p_details():
	
	# default (str): Define what kind of parameters should be used for the simulation or inverse problem ('k': rate constants, 'g': free energy, 'a': activation, 'l': linked)
	
	## Does the lower bound make sense?

	def __init__(self):
		self._k = 0
		self._Ao = None
		self._Ea = None
		self._Ga = None
		self._dG = None
		self._link = None
		self._use = 'k'
		self._lower_bound = None
		self._upper_bound = None

	@property
	def k(self):
		return self._k
	@k.setter
	def k(self,a):
		self._k = a
		self._use = 'k'

	@property
	def Ao(self):
		return self._Ao
	@Ao.setter
	def Ao(self,a):
		self._Ao = a
		self._use = 'E'

	@property
	def Ea(self):
		return self._Ea
	@Ea.setter
	def Ea(self,a):
		self._Ea = a
		self._use = 'E'

	@property
	def Ga(self):
		return self._Ga
	@Ga.setter
	def Ga(self,a):
		self._Ga = a
		self._use = 'G'

	@property
	def dG(self):
		return self._dG
	@dG.setter
	def dG(self,a):
		self._dG = a
		self._use = 'G'

	@property
	def use(self):
		return self._use
	@use.setter
	def use(self,a):
		self._use = a

	@property
	def link(self):
		return self._link
	@link.setter
	def link(self,a):
		self._link = a
		self._use = 'l'

	@property
	def lower_bound(self):
		return self._lower_bound
	@lower_bound.setter
	def lower_bound(self,a):
		self._lower_bound = a

	@property
	def upper_bound(self):
		return self._upper_bound
	@upper_bound.setter
	def upper_bound(self,a):
		self._upper_bound = a