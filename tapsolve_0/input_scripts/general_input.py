from tap_sim import tap_simulation_function
from math import *
import numpy as np
import time
import sys
import matplotlib.pyplot as plt
	
######## Reactor type & data storage ########################################################################

fig, ax = plt.subplots()
ax.set_ylabel('$Dimensionless\ Flux\ (1/s)$')
ax.set_xlabel('$Time\ (s)$')

colors = ['b','r','g','m','k']

reac1 = ['A + * -> A*']

reac2 = ['A + * <-> A*']

reac3 = ['A + * -> A*','B + * -> B*']

reac4 = ['A + * <-> A*','B + * <-> B*'] 

reac5 = ['A + * <-> A*','B + * <-> B*','C + * <-> C*','A* + B* <-> C*'] 

reac6 = ['A + * <-> A*','B + * <-> B*','C + * <-> C*','D + * -> D*','A* + B* <-> C*','A* + C* -> D*' ]
#
reactions_list = [reac1,reac2,reac3,reac4,reac5,reac6]
#
#mesh_size_list = [100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500]
time_step_list = [2,20,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000]

store_data_csv = './variation_data.csv'

#output_data_simulation = np.asarray(mesh_size_list)
output_data_simulation = np.asarray(time_step_list)

reactants_num_list = [1,1,2,2,2,2]
reactant_ratio_list = [[1],[1],[1,1],[1,1],[1,1],[1,1]]
rangesurface_species_list = [[0,(0.5/1000)*6.022e23],[0,(0.5/1000)*6.022e23], [0,0,(0.5/1000)*6.022e23], [0,0,(0.5/1000)*6.022e23], [0,0,0,(0.5/1000)*6.022e23],[0,0,0,0,(0.5/1000)*6.022e23] ]

for what_a_concept_num,what_a_concept in enumerate(reactions_list):
	sim_time_list = []
	for k_cool in range(0,len(time_step_list)):	
	#for time_step_value in time_step_list:	
		time_step_value = 600
		print("CURRENT VALUE: "+str(time_step_list[k_cool]))
		reactor_type = 'tap' 	#currently t_pfr / tap / t_pfr_diff / batch {not ready just yet}
		output_file_name = "test" # "./FILENAME.csv"
		theta = 1 				# forward_euler = 0, backward_euler = 1, crank_nicolson = 1/2
		solver_method_options = 'None'	#'simple_adaptive' or None 
		Time = 0.4				# Total pulse time
		Time_steps = time_step_list[k_cool]  #time_step_value	# Number of time steps  
		mesh_size = 280#mesh_size_list[k_cool]#mesh_size_list[k_cool] #280	# Want to test and store all of these
		T = 400 
		sensitivity_analysis = False # True = on / False = off
		Frequency_of_sensitivity = 1
		Display_figure = False
		save_figure = False
		store_data = False
		
		############# Feed composition ############################################################################
		
		reactants_num = reactants_num_list[what_a_concept_num]
		Inert_pulse_size = (5e-8)*6.022e23	###Size of inert pulse (# of molecules, not mol)#5e-8 #(5e-9)*
		reactant_ratio = reactant_ratio_list[what_a_concept_num]		###Ratio of reactant to the size of Inert_pulse_size
		number_of_pulses = 1#reactants_num_list[what_a_concept_num]
		
		############# Reaction Equations ############################################################################
		
		reactions_test = what_a_concept

		### Pt oxidation
		# actual
		#reactions_test = ['O2 + 2* -> 2O*']
		# Toy model
		#reactions_test = ['A + * -> A*']
		
		### CO OXIDATION
		#Eley
		#reactions_test = ['A + * -> A*','B + A* -> * + C']
		#Eley with preadsorbed
		#reactions_test = ['CO + O* ->  CO2 + *']
		#Lang
		#reactions_test = ['CO + * <-> CO*','CO* + O* -> CO2 + 2*']
		#Eley & Lang
		#reactions_test = ['CO + * <-> CO*','CO* + O* -> CO2 + 2*','CO + O* -> CO2 + *']
		#El_La_aop
		#reactions_test = ['CO + * <-> CO*','CO* + O* -> CO2 + 2*','CO + O* -> CO2 + *','O* + z* -> zO* + *']
		#Reece reaction network
		#reactions_test = ['CO <-> CO*', 'CO* + OA* -> CO2*', 'CO* + OB* -> CO2*', 'CO2* -> CO2']
		
		
		### Ammonia decomposition 
		#reactions_test = ['NH3 + * <-> NH3*', 'H2 + 2* <-> 2H*', 'N2 + * <-> N2*', 'N2* + * <-> 2N*', 'NH3* + * <-> NH2* + H*', 'NH2* + * <-> NH* + H*', 'NH* + * <-> N* + H*']
		
		############# Reactor Parameters ############################################################################
		
		len_inert_1 = 2.7/2#2.54*(19/20)/2		###Length of the first inert zone (cm)
		len_cat = 0.1#2.54*(1/20)			###Length of the catalyst zone (cm) 
		len_inert_2 = 2.7/2#2.54*(19/20)/2		###Length of the second inert zone (cm)
		reac_radius = sqrt((1.8e-5)*(100**2)/3.14159)#sqrt((1.8e-5)*(100**2)/3.14159)		###Radius of the reactor
		
		############# Transport properties ##########################################################################
		
		void_inert = 0.53		###Voidage of the catalyst
		void_cat = 0.53			###Voidage of the catalyst
		
		ref_rate_inert = 43.5319		###Reference 
		ref_rate_cat = 32.5319
		ref_mass = 423.15
		#mass_list = np.array((16,16))
		mass_list = np.array((16,16,16,16,16))
		#mass_list = np.array((28,44,40))	
		#mass_list = np.array((17.03,2.015,28.01,40))
		velocity = 1
		
		############# Catalyst properties ##########################################################################
		
		#Fraction of sites
		mcat = 0.025   			###Catalyst mass in grams
		ssarea = 5E4			###Surface area of support in cm2/g
		bsarea = ssarea*mcat	###Total surface area of support
		sdensity = 1.4e15	###Surface Density (atoms/cm2)
		atloading = 0.08		### Fractional metal loading of catalyst (assume same surface density as support)
		sarea = bsarea*atloading	### Total surface area of metal
		open_sites = 0.55*sdensity
		OA = 0.15*sdensity 		### OA site density (atoms/cm2)
		OB = 0.30*sdensity	
		
		rangesurface_species = [(0.25/3000)*6.022e23,(0.5/1000)*6.022e23]
		#rangesurface_species = [(5e-9)*6.022e23,(5e-9)*6.022e23]
		#rangesurface_species = [((1e7)*(5e-9)/(100**2))*6.022e23]
		#rangesurface_species = [0,OA,0,OB,0]
		
		########## Input Rate Constants #############################################################################
		
		###### Time stepping/solver options
		
		Ke0 = (100**3)/(6.022e23)
		Kd0 = 10
		Ke1 = 400*(1000)/(6.022e23)
		Kd1 = 400*(1000)/(6.022e23)
		Ke2 = (100**3)/(6.022e23)
		Kd2 = 10
		Ke3 = 400*(1000)/(6.022e23)
		Kd3 = 400*(1000)/(6.022e23)
		Ke4 = 400*(1000)/(6.022e23)
		Kd4 = 400*(1000)/(6.022e23)
		Ke5 = 400*(1000)/(6.022e23)
		Kd5 = 400*(1000)/(6.022e23)
		
	
		Ke0 = 1.15*(100**3)/(6.022e23)
		Kd0 = 10
		
		Ke1 = 400*(1000)/(6.022e23)		###Rate constants
		
		Ke2 = (100**3)/(6.022e23)
		Ke3 = 100*(1000)/(6.022e23)
		
		#Ke0 = 1.1108e07
		#Kd0 = 1.6900e12
		#Ke1 = 3.354e-7
		#Ke2 = 3.287e-10
		#Ke3 = 9.6077e9
		
		reactor_kinetics_input = dict(locals())
		
		time_of_sim, graph_data, legend_label = tap_simulation_function(reactor_kinetics_input)
		
		sim_time_list.append(time_of_sim)
	
		for k,j in enumerate(graph_data):
			if j != 'timing':
				if time_step_list[k_cool] > 2:
				#if mesh_size_list[k_cool] > 100:
					ax.plot(graph_data['timing'],graph_data[j],color=colors[k], ls = '--', alpha=0.7)
				else:
					ax.plot(graph_data['timing'],graph_data[j],color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)
				#new_data = np.asarray(graph_data[j])
				#dictionaray_of_numpy_data[legend_label[k-1]] = np.vstack((dictionaray_of_numpy_data[legend_label[k-1]],new_data))
			else:
				pass
		#plt.show()
		reactor_kinetics_input = {}
	temp_data = np.asarray(sim_time_list)
	output_data_simulation = np.vstack((output_data_simulation,temp_data))

np.savetxt(store_data_csv, np.transpose(output_data_simulation), delimiter=",")

ax.legend(title="Gas Species")
plt.show()


####
#eley
# time step
# [0.03116893768310547, 0.15613675117492676, 0.35606884956359863, 0.6849236488342285, 1.0598032474517822, 1.3329894542694092, 1.6916074752807617, 2.0225236415863037, 2.3039438724517822, 2.6381595134735107, 2.977360725402832, 3.294581651687622, 3.5979976654052734, 3.9508378505706787, 4.225225448608398, 4.5638747215271, 4.9195544719696045, 5.208940029144287, 5.486323833465576, 5.858653545379639, 6.077473163604736, 6.385039806365967]
# mesh size
# [2.4852454662323, 2.571659564971924, 2.7827632427215576, 2.8785579204559326, 3.103823184967041, 3.207014322280884, 3.3449018001556396, 3.3801722526550293, 3.408679723739624, 3.518083333969116, 3.6205976009368896, 3.704531192779541, 3.8096961975097656, 3.9312050342559814, 4.061707019805908, 4.1477437019348145, 4.411125421524048, 4.6940789222717285, 5.325071811676025, 5.107831239700317, 5.2876410484313965, 5.328880786895752, 5.785737037658691, 5.151878833770752, 5.327051401138306]

#lang
# time step
# [0.03689098358154297, 0.19363927841186523, 0.4724266529083252, 0.9303324222564697, 1.3787055015563965, 1.7528800964355469, 2.1484787464141846, 2.601935625076294, 2.962425708770752, 3.3792712688446045, 3.720546007156372, 4.056863307952881, 4.504999399185181, 4.926687717437744, 5.293976068496704, 5.750283241271973, 6.869592189788818, 7.668954849243164, 7.488609790802002, 7.574548244476318, 7.805788278579712, 8.134341478347778]
# mesh size
# [2.692828416824341, 2.782825231552124, 3.1577277183532715, 3.318420171737671, 3.4715986251831055, 3.6013758182525635, 3.7780892848968506, 3.9362659454345703, 4.049852609634399, 4.224122762680054, 4.427903413772583, 4.530793905258179, 4.726433753967285, 4.895007848739624, 4.9808030128479, 5.196974039077759, 5.432405710220337, 5.430123567581177, 5.5697925090789795, 5.731900453567505, 6.298583030700684, 6.324296951293945, 6.428957939147949, 6.604550361633301, 6.761291265487671]

#lang + eley
# time step
# [0.04178595542907715, 0.19973134994506836, 0.47632479667663574, 0.9451272487640381, 1.3795580863952637, 1.770157814025879, 2.1658120155334473, 2.6959831714630127, 2.966167449951172, 3.437129497528076, 3.7810492515563965, 4.177079916000366, 4.6083385944366455, 5.020092248916626, 5.415220499038696, 5.82276177406311, 6.169091701507568, 6.600710153579712, 6.913733005523682, 7.348866701126099, 7.768005609512329, 8.111349821090698]
# mesh size
# [2.778442144393921, 2.883129596710205, 3.2374422550201416, 3.399101734161377, 3.5660641193389893, 3.6870150566101074, 3.86112117767334, 4.002504587173462, 4.130796670913696, 4.296316385269165, 4.469245433807373, 4.563286542892456, 4.721400737762451, 4.92127537727356, 5.255516052246094, 5.39132833480835, 5.533004522323608, 5.5011146068573, 5.7149293422698975, 5.994065761566162, 6.121885061264038, 6.148864269256592, 6.377568244934082, 6.443218469619751, 6.757220506668091]

# Reece_CO_oxidation
# time step
# [0.13676238059997559, 0.4782297611236572, 1.1897130012512207, 2.0325169563293457, 2.988412380218506, 3.6358511447906494, 4.695979595184326, 5.597567558288574, 6.403910160064697, 7.399029970169067, 8.061829328536987, 8.935603618621826, 10.205981969833374, 11.194551944732666, 11.59661316871643, 12.542011976242065, 13.925042152404785, 14.802717924118042, 16.15423846244812, 16.871363639831543, 17.149802684783936, 18.342536211013794]
# mesh size
# [5.220643520355225, 5.625455379486084, 6.2563230991363525, 6.398206949234009, 6.7022857666015625, 7.037503719329834, 7.136411428451538, 9.108701705932617, 9.981152534484863, 10.11208724975586, 10.91795563697815, 10.887830257415771, 11.65386414527893, 12.085758924484253, 12.542765617370605, 14.214283227920532, 13.914568424224854, 14.912169218063354, 17.001413822174072, 17.8608295917511, 15.657803773880005, 16.618581771850586, 17.9334557056427, 18.816773891448975, 19.04019260406494]
