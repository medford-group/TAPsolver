from tapsolver import *
import os
import scipy.stats as sp

opdh_2 = read_TAPobject('./mech_1_run.json')

print(opdh_2.output_name)
print('divide')
forward_problem(1.5,10,opdh_2)
flux_graph(opdh_2)
temp_data = read_experimental_data_object(opdh_2.output_name+'/TAP_experimental_data.json')

print(temp_data['C3H8'][1][0])
print(temp_data['C3H8'][2][0])
print(temp_data['C3H8'][3][0])

new_value = sp.moment(temp_data['C3H8'][0],moment=0)
print(new_value)
new_value = sp.moment(temp_data['C3H8'][1],moment=0)
print(new_value)
new_value = sp.moment(temp_data['C3H8'][2],moment=0)
print(new_value)
new_value = sp.moment(temp_data['C3H8'][3],moment=0)
print(new_value)
new_value = sp.moment(temp_data['C3H8'][4],moment=0)
print(new_value)

#new_value = [sp.moments(moment=0),sp.moments(moment=1),sp.moments(moment=2)]
flux_graph(opdh_2)
#os.rmdir('./exp_new')
