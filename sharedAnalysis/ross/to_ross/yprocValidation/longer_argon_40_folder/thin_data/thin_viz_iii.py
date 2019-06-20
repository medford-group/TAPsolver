import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp
import imageio
import sys
import os

sim_info = pd.read_csv('../input_file.csv',header = None)
curr_time_steps = 1000

length = float(sim_info[1][3])
mesh_size = float(sim_info[1][4])
time_steps = float(sim_info[1][2])
time_tot = float(sim_info[1][1])
cat_frac = float(sim_info[1][5])
reac_radius = float(sim_info[1][6])

kf = 10
kb = 30

gas_data = pd.read_csv('./'+'CO'+'.csv',header = None)
#print(thin_data.iloc[2].sum()*)
surf_data = pd.read_csv('./'+'CO*'+'.csv',header = None)
open_data = pd.read_csv('./'+'*'+'.csv',header = None)

x_length = mp.ceil(cat_frac*mesh_size)

x_dim = []




for k in range(x_length+1):
	x_dim.append((k+1)*length/mesh_size)


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


fig, host = plt.subplots()
fig.subplots_adjust(right=0.75)

par1 = host.twinx()
par2 = host.twinx()

par2.spines["right"].set_position(("axes", 1.2))

make_patch_spines_invisible(par2)

par2.spines["right"].set_visible(True)


for j in range(2,7):
	time = []
	av_conc = []
	av_surf = []
	av_rate = []
	for step in range(curr_time_steps):
		av_cv = gas_data.iloc[step][j]#.sum()*(length/mesh_size)/(length*cat_frac)
		av_conc.append(av_cv)
		av_sv = surf_data.iloc[step][j]#.sum()*(length/mesh_size)/(length*cat_frac)
		av_surf.append(av_sv*(9.6*10**(-4)))
		av_op = open_data.iloc[step][j]#.sum()*(length/mesh_size)/(length*cat_frac)
		av_rate.append( (-kb*av_sv + kf*av_cv*av_op)*(9.6*10**(-4)) )
		time.append(step*(time_tot/time_steps))

	p1, = host.plot(time, av_rate, "b-")
	p2, = par1.plot(time,av_conc, "g-")
	p3, = par2.plot(time, av_surf, "r-")


host.set_xlim(-0.02, 0.52)
host.set_ylim(-.3, 4)
par1.set_ylim(-.3, 4)
par2.set_ylim(-0.015, 0.2)


host.set_xlabel("$t\ (s)$")
host.set_ylabel("$R\ (mmol/kg/s)$")
par1.set_ylabel("$Cg\ (mmol/m3)$")
par2.set_ylabel("$Cz\ (mmol/kg)$")



host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
par2.yaxis.label.set_color(p3.get_color())


tkw = dict(size=4, width=1.5)
host.tick_params(axis='y', colors=p1.get_color(), **tkw)
par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
host.tick_params(axis='x', **tkw)

host.hlines(y=0,xmin=0,xmax=0.5,linestyle='--',colors = 'k')


rateData = pd.DataFrame(av_rate) 
rateData.to_csv('./co_rate.csv',header=None)		

plt.savefig('./me_figure.png')

plt.show()

