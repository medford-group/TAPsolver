#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 12:54:23 2023

@author: ayonge
"""

from .p_details import p_details
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

class process(object):

	"""
	
	This class provides all kinetic information for an elementary process.
	
	Args:
		
		elementary_processes (int): The number of cells in the mesh used to solve the PDEs.

		rate_array (numpy array): The stoichiometric matrix of the microkinetic model.

		reactions (list): The individual elementary processes introduced in the mechanism (also found in elementary_processes).

		reactants (list): The reactive gas and surface species in the mechanism.   

	"""

	def __init__(self,processString):
		self.reaction = processString
		self.f = p_details()
		self.b = p_details()