import numpy as np
from scipy.optimize import minimize

z_sl_LA = np.array([0.19, 0.13, 0.86])
z_sl_BA = np.array([0.24, 0.34, 0.40]) 

PO_sl_LA = np.array([0.04, 0.46, 0.28])
PO_sl_BA = np.array([0.07, 0.41, 0.4])

sl_LA = z_sl_LA
sl_BA = z_sl_BA
S = 1.33 #cm2
def objective(e, LA, BA, surface):
	#print(e)
	diff = np.empty([sl_LA.size])
	for j in range(len(sl_LA)):
		tempValue = - S * ((sl_LA[j]/e[0])+(sl_BA[j]/e[1]))
		if abs(tempValue) < 1:
			diff[j] = 1 + tempValue
		else:
			diff[j] = 0
	resid = np.sum(diff**2)
	print(resid)
	return(resid)
	
Emeis = [1.67, 2.22]
Datka = [1.11, 0.73]
Tamura = [1.23, 1.73]

e_init = Emeis
objective(e_init, sl_LA, sl_BA, S)

res = minimize(objective, e_init, args=(sl_LA, sl_BA, S), tol=1e-20)
print(res)

print(res.x)
print('initial epsilon = ', e_init)
print('final epsilon = ', res.x)