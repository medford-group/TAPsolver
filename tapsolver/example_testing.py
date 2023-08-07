from data import *
#from experimental_data import experimental_data
from .mechanism import *#mechanism
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

no.reactor.radius = 0.2

no.reactor.lengths = {0:1.85,1:0.36,2:1.85}

no.reactor.reference_diffusions = {0:16,1:16,2:16}

#no.reactor.lengths = {0:1,1:0.1,2:1,3:0.1,4:1}

no.mechanism.processes = [0,'A + *(0) <-> A*(0)']
no.mechanism.processes = [1,'A*(0) <-> B*(0)']
no.mechanism.processes = [2,'A + B*(0) <-> C*(0)']
no.mechanism.processes = [3,'B*(0) <-> B + *(0)']
no.mechanism.processes = [4,'C*(0) <-> C + *(0)']
no.mechanism.processes = [5,'C*(0) <-> D*(0)']

no.mechanism.processes[0].f.k = 80#5
no.mechanism.processes[1].f.k = 7#100
no.mechanism.processes[2].f.k = 50#100
no.mechanism.processes[3].f.k = 10#5
no.mechanism.processes[4].f.k = 14#100
no.mechanism.processes[5].f.k = 1#100
## Automate this so user doesn't have to call thisred
no.species.initialize(no.mechanism.reactants)

no.species.gasses['A'].mass = 24
no.species.gasses['B'].mass = 19
no.species.gasses['C'].mass = 14
no.species.gasses['A'].noise = 0.0
no.species.gasses['A'].intensity = 1

#no.species.gasses['B'].noise = 0.0

argon = gas(40)
argon.intensity = 1
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

no.species.adspecies['*(0)'].conc = 10
no.species.adspecies['*(0)'].zones = [1]

no.display_analytical = False

no.output_name = 'test_fit_data'
##no.data_name = 'test_exp_data/TAP_experimental_data.json'

no.refine = {0:{0:[no.reactor.fractions[0],no.reactor.fractions[0]+no.reactor.fractions[1]],1:2}}

no.gasses_objective = ['A','B','C']
no.poi = ['[0].f.k','[1].f.k','[2].f.k','[3].f.k','[4].f.k']

no.optimize = False

#no.view_mesh = True

no.conc_profiles = ['A','C','D*(0)','*(0)']
#no.conc_profiles = ['D*(0)']

#forward_problem(1.5,50,no)

#moment_graph(no)

#generateGif(no)
##fluxGif(no)
single_pulse_concentrationGif(no,0)
#concentrationGif(no,43)
#sys.exit()
#forward_problem(1.5,50,no)

#no.pulses_graphed = range(10,20)

##flux_graph(no)

moment_graph(no)