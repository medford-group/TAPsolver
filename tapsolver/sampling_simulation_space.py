import scipy.stats as sp
# load tapobject
read_TAPobject('./quick_test.json')
# run simulation
forward_problem(5,1,opdh_1)
# load data
read_experimental_data_object(file_name)
# calculate moments
new_value = [sp.moments(moment=0),sp.moments(moment=1),sp.moments(moment=2)]
# store moments 

# delete simulated data