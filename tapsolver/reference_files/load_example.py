
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pickle
import os

def load_example(reactor_name):
	
	this_dir, this_filename = os.path.split(__file__)
	new_reactor = pickle.load(open(this_dir+'/'+reactor_name,'rb'))
	
	return new_reactor
	