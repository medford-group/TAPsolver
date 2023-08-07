from reactor import *
from species import *
from mechanism import *

from TAPobject import *
from forward_problem import *
import sys

nt = TAPobject()

print(nt.reactor.lengths)

nt.mechanism.processes[0] = process("A + *(0) <-> A*(0)")
nt.mechanism.processes[1] = process("B + *(0) <-> B*(0)")
nt.mechanism.processes[2] = process("A*(0) + B*(0) <-> C + 2*(0)")

mechanism_constructor(nt.mechanism)

nt.species.initialize(nt.mechanism.reactants)

nt.species.gasses["A"].intensity = 2
nt.species.gasses["A"].mass = 32
nt.species.gasses["B"].intensity = 1.5
nt.species.gasses["B"].mass = 28
nt.species.gasses["C"].intensity = 0
nt.species.gasses["C"].mass = 60

nt.mechanism.write_mech_excel('test_store',rewrite=False,use_data=True)
nt.mechanism.store_matrix('test_store_matrix')

forward_problem(1,1,nt)

sys.exit()

nr = reactor()

nm = mechanism()

nm.add(0,process("CH3OCH3 + H*(0) <-> CH3OCH3H*(0)"))

nm.processes[0].f.k = 0
nm.processes[0].b.k = 0

#nm.processes = temp

#mechanism_constructor(nm)

nm.add(1,process("CH3OCH3 + H*(1) <-> CH3OCH3H*(1)"))

