import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
fig,ax = plt.subplots()

a = np.random.random((4, 4))
for j1 in range(0,4):
	for j2 in range(0,4):
		with open('./'+str(j1)+'_'+str(j2)+'_'+str(0)+'_'+str(0)+'/exp_new/optimal_criteria.txt') as f:
			lines = f.readlines()
		f.close()
		a[j1][j2] = lines[1]


ax = sns.heatmap(a, linewidth=0.5)
plt.xlabel('test name')
plt.ylabel('test name')
#plt.imshow(a, cmap='hot', interpolation='nearest')
plt.show()