from fenics import *
from fenics_adjoint import *
from structures import *
from reactor_species import *
from file_io import *
from mechanism_construction import *
from forward_problem import *
import jsonpickle
import json
import sys
import copy

# Add data storage
# Add thin zone storage
# Add data read
# Add object storage
# Add objective function reading
# Add control specifications
# Read / make new objects

# 1. Edit TAP object
# 2. Store in the new folder
# 3. Edit TAP object

testGen1 = reactor()
testGen1.zone_lengths = {0: 2.5, 1: 0.1, 2: 2.5}
testGen1.zone_voids = {0: 0.4, 1: 0.4, 2: 0.4}
testGen1.reactor_radius = 1.0

#testGen1 = readCSV_reactor('./input_file_2.csv')
new_reactor_species = reactor_species()

new_reactor_species.inert_diffusion = 16
new_reactor_species.catalyst_diffusion = 16
new_reactor_species.reference_temperature = 385.6
new_reactor_species.reference_mass = 40
new_reactor_species.temperature = 873
new_reactor_species.reference_pulse_size = 1
new_reactor_species.advection = 0

C3H8 = define_gas()
C3H8.mass = 44
C3H8.intensity = 1
C3H8.delay = 0.0
new_reactor_species.add_gas('C3H8',C3H8)

C3H6 = define_gas()
C3H6.mass = 42
C3H6.intensity = 0
C3H6.delay = 0.0
new_reactor_species.add_gas('C3H6',C3H6)

H2 = define_gas()
H2.mass = 2
H2.intensity = 0
H2.delay = 0.0
new_reactor_species.add_gas('H2',H2)

argon = define_gas()
argon.mass = 44
argon.intensity = 1
argon.delay = 0.0
new_reactor_species.add_inert_gas('argon',argon)

#jsonStr = json.dumps(vars(new_reactor_species))

#sys.exit()

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H8*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H7*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('H^',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H6*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('H2^',s)

s = define_adspecies()
s.concentration = 0.7
new_reactor_species.add_adspecies('*',s)

s = define_adspecies()
s.concentration = 0.7
new_reactor_species.add_adspecies('^',s)

new_mechanism = mechanism()

new_mechanism.elementary_processes[0] = elementary_process()
new_mechanism.elementary_processes[0].processString = 'C3H8 + * <-> C3H8*'
new_mechanism.elementary_processes[0].forward = elementary_process_details()
new_mechanism.elementary_processes[0].backward = elementary_process_details()

new_mechanism.elementary_processes[1] = elementary_process()
new_mechanism.elementary_processes[1].processString = 'C3H8* + ^ <-> C3H7* + H^'
new_mechanism.elementary_processes[1].forward = elementary_process_details()
new_mechanism.elementary_processes[1].backward = elementary_process_details()

new_mechanism.elementary_processes[2] = elementary_process()			
new_mechanism.elementary_processes[2].processString = 'C3H7* + ^ <-> C3H6* + H^'
new_mechanism.elementary_processes[2].forward = elementary_process_details()
new_mechanism.elementary_processes[2].backward = elementary_process_details()

new_mechanism.elementary_processes[3] = elementary_process()			
new_mechanism.elementary_processes[3].processString = 'C3H6* <-> C3H6 + *'
new_mechanism.elementary_processes[3].forward = elementary_process_details()
new_mechanism.elementary_processes[3].backward = elementary_process_details()

new_mechanism.elementary_processes[4] = elementary_process()			
new_mechanism.elementary_processes[4].processString = '2H^ <-> H2^'
new_mechanism.elementary_processes[4].forward = elementary_process_details()
new_mechanism.elementary_processes[4].backward = elementary_process_details()

new_mechanism.elementary_processes[5] = elementary_process()			
new_mechanism.elementary_processes[5].processString = 'H2^ <-> H2 + 2^'
new_mechanism.elementary_processes[5].forward = elementary_process_details()
new_mechanism.elementary_processes[5].backward = elementary_process_details()

mechanism_constructor(new_mechanism)

new_mechanism.elementary_processes[0].forward.k = 1.66e8
new_mechanism.elementary_processes[0].backward.k = 8.23e15
new_mechanism.elementary_processes[1].forward.k = 9.71e11
new_mechanism.elementary_processes[1].backward.k = 6.63e12
new_mechanism.elementary_processes[2].forward.k = 2.25e12
new_mechanism.elementary_processes[2].backward.k = 1.7e11
new_mechanism.elementary_processes[3].forward.k = 1.02e15
new_mechanism.elementary_processes[3].backward.k = 1.69e8
new_mechanism.elementary_processes[4].forward.k = 7.25e12
new_mechanism.elementary_processes[4].backward.k = 5.16e10
new_mechanism.elementary_processes[5].forward.k = 1.23e13
new_mechanism.elementary_processes[5].backward.k = 5.64e8

TAP_test = TAPobject()
TAP_test.mechanism = new_mechanism
TAP_test.reactor_species = new_reactor_species
TAP_test.reactor = testGen1

forward_problem(1,20,TAP_test)
