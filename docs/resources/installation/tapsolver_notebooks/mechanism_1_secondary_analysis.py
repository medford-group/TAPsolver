from tapsolver import *

mech_1 = read_TAPobject('./mech_1.json')

mech_1_up = update_parameters(mech_1)

mech_1_up.uncertainty = True

forward_problem(0.3,1,mech_1_up)