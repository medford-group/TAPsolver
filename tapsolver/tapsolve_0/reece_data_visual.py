import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

temp = pd.read_csv("./reece_output_data.csv")
time = temp.ix[:,0]
inert = temp.ix[:,1]
CO = temp.ix[:,2]
CO2 = temp.ix[:,3]
#experimental_data[new_name] = temp.ix[:,desired_column]
#temp = pd.read_csv(file_location+i+"ReactionRate.csv")
#y_column['y'] = temp.ix[:,desired_column]

fig2, ax2 = plt.subplots()
ax2.set_xlabel('$t (s)$')
ax2.set_ylabel('$Flux (molecules/s)$')

legend = ['CO','inert','CO2']
colors = ['b','r','m','g','o']
ax2.plot(time,inert,color=colors[0],label=legend[0],ls = '--', alpha=0.7)
ax2.plot(time,CO,color=colors[1],label=legend[1],ls = '--', alpha=0.7)
ax2.plot(time,CO2,color=colors[2],label=legend[2],ls = '--', alpha=0.7)

ax2.legend(title="Gas Species")
plt.show()
