import numpy as np
​
z_sl_LA = np.array([0.19, 0.13, 0.86]) #umol cm-1 - H-beta (25) H-ZSM-5 (30) H-beta (38) - (sio2/al2o3)  - 
z_sl_BA = np.array([0.24, 0.34, 0.40]) 
​
# PHOSPHATES
PO_sl_LA = np.array([0.04, 0.46, 0.28]) #umol cm-1 - ZrPO NbPO LaPO
PO_sl_BA = np.array([0.07, 0.41, 0.4])
​
sl_LA = z_sl_LA
sl_BA = z_sl_BA
S = 1.33 #cm2
print(z_sl_LA.size)
def objective(e, LA, BA, surface):
	print(e)
	diff = np.empty([sl_LA.size])
	for j in range(len(sl_LA)):
		diff[j] = 1 - S * (sl_LA[j]/e[0]+sl_BA[j]/e[1])
	print(diff)
	resid = np.sum(diff**2)
	print(resid)
	return(resid)
Emeis = [1.67, 2.22] # [LA, BA]
Datka = [1.11, 0.73] # [LA, BA]
Tamura = [1.23, 1.73] # [LA, BA]
from scipy.optimize import minimize
e_init = Emeis
objective(e_init, sl_LA, sl_BA, S)
​
eBOUNDS = ( (0,np.inf), (0,np.inf) )
res = minimize(objective, e_init, args=(sl_LA, sl_BA, S), bounds=eBOUNDS, tol=1e-20)
print(res)
​
print(res.x)
print('initial epsilon = ', e_init)
print('final epsilon = ', res.x)