import math as mp

kbt = 1.38064852e-23
hb = 6.62607004e-34
Rgas = 8.314
T = 340


test = 1
#test = (kbt*T/hb)*mp.exp(-G/(Rgas*T))

G = -mp.log(test*hb/(kbt*T))*Rgas*T

G1f = 83640.78425382517/2826.76
G1b = 90149.63969129704/2826.76
G2f = 90149.63969129704/2826.76
G2b = 91157.87415582528/2826.76
G3f = 70623.07337888148/2826.76
G3b = 109676.20600371258/2826.76

print(G1f)
print(G2f)
print(G3f)

print(G1f - G1b)
print(G2f - G2b)
print(G3f - G3b)

test1 = (kbt*T/hb)*mp.exp(-G1f/(Rgas*T))
test2 = (kbt*T/hb)*mp.exp(-G2f/(Rgas*T))
test3 = (kbt*T/hb)*mp.exp(-G3f/(Rgas*T))

print(test1)
print(test2)
print(test3)

print((kbt*T/hb)*mp.exp((-G1f)/(Rgas*T)))
print((kbt*T/hb)*mp.exp((-G2f)/(Rgas*T)))
print((kbt*T/hb)*mp.exp((-G3f)/(Rgas*T)))