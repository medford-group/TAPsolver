
from .mechanism import *
import numpy as np

def kinetics_sheet(mechanism_data: mechanism):
	print(len(mechanism_data.reactions))
	print(mechanism_data.processes[0].f.__dict__)