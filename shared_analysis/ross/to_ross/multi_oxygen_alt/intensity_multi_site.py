import pandas as pd
import matplotlib.pyplot as plt
import sys

fig2,ax2 = plt.subplots()

peak = []
value = []

user_data = pd.read_csv('./lang_eley/e_n3_folder/flux_data/CO.csv',header=None)

peak_loc = user_data[len(user_data.columns)-1].iloc[user_data[len(user_data.columns)-1].idxmax()]

final_peak = peak_loc

for k in range(0,len(user_data.columns)-1):
	peak.append(k+1)
	peak_loc = user_data[k+1].iloc[user_data[k+1].idxmax()]
	#print(peak_loc)
	value.append((final_peak - peak_loc)/final_peak)

ax2.scatter(peak,value)

ax2.set_xlabel('Pulse Number')

ax2.set_ylabel('Relative CO Pulse Intensity')
plt.savefig('./lang_eley/e_n3_folder/multi_site_example.png')
plt.show()




