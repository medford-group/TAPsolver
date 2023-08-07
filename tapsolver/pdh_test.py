from data import *
#from experimental_data import experimental_data
from mechanism import *#mechanism
from reactor import *
from species import *
from TAPobject import *#TAPobject
from read_files import *
from etc import *
from flux_graph import *
from moment_graph import *

import sys

from forward_problem import *

no = TAPobject()

no.reactor.radius = 0.19

no.reactor.lengths = {0:1.85,1:0.36,2:1.85}

print(no.reactor.reference_diffusions)
no.reactor.reference_temperature = 700
no.reactor.reference_mass = 40
no.reactor.voids = {0:0.4,1:0.4,2:0.4}
no.species.temperature = 833
no.reactor.temperature = 833
no.reactor.reference_diffusions = {0:25.04,1:25.04,2:25.04}
print(no.reactor.reference_diffusions)

#no.reactor.lengths = {0:1,1:0.1,2:1,3:0.1,4:1}

no.mechanism.processes = [0,'C3H8 + *(0) <-> C3H8*(0)']
no.mechanism.processes = [1,'C3H8*(0) + *(0) <-> C3H6*(0) + H2*(0)']
no.mechanism.processes = [2,'C3H6*(0) <-> C3H6 + *(0)']
no.mechanism.processes = [3,'H2*(0) <-> H2 + *(0)']
no.mechanism.processes = [4,'C3H8*(0) + *(0) <-> C2H5*(0) + CH3*(0)']
no.mechanism.processes = [5,'C3H6*(0) + *(0) <-> C2H5*(0) + CH2*(0)']
no.mechanism.processes = [6,'H2*(0) <-> 2H*(0)']


#no.mechanism.processes[0].f.k = 1.429e-1#5
#no.mechanism.processes[0].b.k = 2.215e-7#5
#no.mechanism.processes[1].f.k = 1.31#100
#no.mechanism.processes[1].b.k = 0.0#100
#no.mechanism.processes[2].f.k = 3.57e2#50
#no.mechanism.processes[2].b.k = 2.46e-5#50
#no.mechanism.processes[3].f.k = 4.20e2#30
#no.mechanism.processes[3].b.k = 1.36e-5#30
#no.mechanism.processes[4].f.k = 5.49e-1#10
#no.mechanism.processes[4].b.k = 0#10
#no.mechanism.processes[5].f.k = 2.96e-1#1
#no.mechanism.processes[5].b.k = 0#1
#no.mechanism.processes[6].f.k = 2.40e-5#1
#no.mechanism.processes[6].b.k = 1.48#1

no.mechanism.processes[0].f.k = 1.429e-5#5
no.mechanism.processes[0].b.k = 2.215e-7#5
no.mechanism.processes[1].f.k = 1.31#100
no.mechanism.processes[1].b.k = 0.0#100
no.mechanism.processes[2].f.k = 3.57e2#50
no.mechanism.processes[2].b.k = 2.46e-5#50
no.mechanism.processes[3].f.k = 4.20e2#30
no.mechanism.processes[3].b.k = 1.36e-5#30
no.mechanism.processes[4].f.k = 5.49e-8#10
no.mechanism.processes[4].b.k = 0#10
no.mechanism.processes[5].f.k = 2.96e-6#1
no.mechanism.processes[5].b.k = 0#1
no.mechanism.processes[6].f.k = 2.40e-5#1
no.mechanism.processes[6].b.k = 1.48#1

## Automate this so user doesn't have to call thisred
no.species.initialize(no.mechanism.reactants)

no.species.gasses['C3H8'].mass = 44.1
no.species.gasses['C3H8'].intensity = 6.87

no.species.gasses['C3H6'].mass = 42.08

no.species.gasses['H2'].mass = 2

argon = gas(40)
argon.intensity = 6.87
no.species.add_inert_gas('Ar',argon)
print(no.species.inert_gasses['Ar'].intensity)
#sys.exit()
#add_inert_gas(self,name='', gas_data = gas)

#no.species.gasses['B'].mass = 8
#no.species.gasses['B'].delay = 0.2
#no.species.gasses['B'].intensity = 0
#
#no.species.gasses['C'].mass = 17
#no.species.gasses['C'].intensity = 0

no.species.adspecies['*(0)'].conc = 2400000
no.species.adspecies['*(0)'].zones = [1]

no.display_analytical = True

no.conc_profiles = ['*(0)']

forward_problem(0.8,10,no)

#no.pulses_graphed = range(10,20)

##flux_graph(no)

moment_graph(no)