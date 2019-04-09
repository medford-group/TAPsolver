import pandas as pd
import matplotlib.pyplot as plt
import sys

df1 = pd.read_csv('./time_sensitivity_data/time_eley_lang.csv',header=None)
df2 = pd.read_csv('./time_sensitivity_data/time-2.csv',header=None)
colors = ['r','b']
legend_label = ['E.R./L.H.','E.R./L.H. w/ A ']

fig2, ax2 = plt.subplots()
ax2.set_xlabel('Time Step')
ax2.set_ylabel('Required Time (min)')
ax2.plot(df1.iloc[:,0],df1.iloc[:,1]/60,color=colors[0],label=legend_label[0], ls = '--', alpha=0.7)
ax2.plot(df2.iloc[:,0],df2.iloc[:,1]/60,color=colors[1],label=legend_label[1], ls = '--', alpha=0.7)
ax2.legend(title="Process")

time_required = []
time_required.append(df1[1].sum()/60/60/24)
time_required.append(df2[1].sum()/60/60/24)

#props = dict(boxstyle='round', facecolor='white', alpha=0.4)

print(time_required)
plt.savefig('./sens_time.png')
plt.show()