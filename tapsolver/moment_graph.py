#from structures import *
from .TAPobject import TAPobject

#from file_io import *
from .read_files import read_experimental_data_object

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
import copy
import sys
import imageio

from matplotlib.backends.backend_agg import FigureCanvasAgg

from .forward_problem import *
def moment_graph(TAPobject_data: TAPobject):

	fig2, ax2 = plt.subplots()

	colors = ['b','orange','g','r','k','y','c','m','brown','darkgreen','goldenrod','lavender','lime']

	def integrate(x, y):
		sm = 0
		for i in range(1, len(x)):
			h = x[i] - x[i-1]
			sm += h * (y[i-1] + y[i]) / 2
		return sm

	
	ax2.set_ylabel('$Moment\ (nmol)$', fontsize = 14)
	ax2.set_xlabel('$Pulse\ Number$', fontsize = 14)

	try:
		#TAPobject_data.experimental_data = #output_data
		experimental_data = read_experimental_data_object(TAPobject_data.data_name)
		experimental_data_exist = True
	except:
		print('no experimental data')
		experimental_data_exist = False
	synthetic_data = read_experimental_data_object('./'+TAPobject_data.output_name+'/TAP_experimental_data.json')
	
	legend_label = []
	for j in TAPobject_data.species.gasses:
		legend_label.append(j)
	for j in TAPobject_data.species.inert_gasses:
		legend_label.append(j)
	gas_moments = {}
	gas_moments['pulse#'] = []
	for jnum,j in enumerate(TAPobject_data.species.gasses):
		gas_moments[j] = []
		if TAPobject_data.pulses_graphed == []:
			for k in synthetic_data[j]:
				if jnum == 0:
					gas_moments['pulse#'].append(int(k))
				gas_moments[j].append(integrate(synthetic_data['time'][0], synthetic_data[j][int(k)]))
				
		else:
			for k in TAPobject_data.pulses_graphed:
				if jnum == 0:
					gas_moments['pulse#'].append(k)
				gas_moments[j].append(integrate(synthetic_data['time'][0], synthetic_data[j][int(k)]))

	for jnum,j in enumerate(TAPobject_data.species.inert_gasses):
		gas_moments[j] = []
		if TAPobject_data.pulses_graphed == []:
			for k in synthetic_data[j]:
				gas_moments[j].append(integrate(synthetic_data['time'][0], synthetic_data[j][int(k)]))
		else:
			for k in TAPobject_data.pulses_graphed:
				gas_moments[j].append(integrate(synthetic_data['time'][0], synthetic_data[j][int(k)]))

	##marker_style = dict(color='tab:blue', linestyle=':', marker='s',markersize=15, markerfacecoloralt='tab:red')

	for jnum,j in enumerate(TAPobject_data.species.gasses):
		plt.scatter(gas_moments['pulse#'],gas_moments[j],label=j,color=colors[jnum],marker='s')

	for jnum,j in enumerate(TAPobject_data.species.inert_gasses):
		plt.scatter(gas_moments['pulse#'],gas_moments[j],label=j,color=colors[len(TAPobject_data.species.gasses.keys())+jnum],marker='s')

	exp_moments = {}
	if experimental_data_exist == True:
		for jnum,j in enumerate(TAPobject_data.species.gasses):#

			if TAPobject_data.pulses_graphed == []:
				
				for k in synthetic_data[j]:
					if jnum == 0:
						gas_moments['pulse#'].append(k)
					exp_moments[j].append(integrate(experimental_data['time'][0][::6], experimental_data[j][int(k)][::6]))
			else:
				for k in TAPobject_data.pulses_graphed:
					if jnum == 0:
						gas_moments['pulse#'].append(k)
					exp_moments[j].append(integrate(experimental_data['time'][0][::6], experimental_data[j][int(k)][::6]))
	
		for jnum,j in enumerate(TAPobject_data.species.inert_gasses):
			for k in experimental_data[j]:
			
				if k == 0:
					plt.scatter(experimental_data['time'][0][::6], experimental_data[j][0][::6],s=5,label=j+'_exp',color=colors[jnum+len(TAPobject_data.species.gasses.keys())])#
				else:
					plt.scatter(experimental_data['time'][0][::6], experimental_data[j][0][::6],s=5,color=colors[jnum+len(TAPobject_data.species.gasses.keys())])
	
	#for jnum,j in enumerate(TAPobject_data.species.gasses):
	#	ax2.plot(exp_moments['pulse#'],exp_moments[j], 'ks', markerfacecolor='none', markeredgecolor=colors[jnum],label='Exp. '+j)

	#for jnum,j in enumerate(TAPobject_data.species.inert_gasses):
	#	ax2.plot(exp_moments['pulse#'],exp_moments[j], 'ks', markerfacecolor='none', markeredgecolor=colors[len(TAPobject_data.species.gasses.keys())+jnum],label='Exp. '+j)
	#	plt.scatter(exp_moments['pulse#'],exp_moments[j],label=j,color=colors[len(TAPobject_data.species.gasses.keys())+jnum],marker='s')

		

	plt.gca().set_ylim(bottom=0)
	plt.savefig('moment_graph.pdf')
	plt.legend()
	plt.show()
	sys.exit()

def generateGif(TAPobject_data: TAPobject):

	"""
	Return a gif showing changes made during the optimization process
	"""
	
	molecules = []
	
	for j in TAPobject_data.species.gasses.keys():
		molecules.append(j)
	
	def add_subplot_axes(ax,rect,axisbg='w'):
		
		"""
		Generates the subplot to help visualize some other quantitiy
		"""

		fig = plt.gcf()
		box = ax.get_position()
		width = box.width
		height = box.height
		inax_position  = ax.transAxes.transform(rect[0:2])
		transFigure = fig.transFigure.inverted()
		infig_position = transFigure.transform(inax_position)    
		
		x = infig_position[0]
		y = infig_position[1]
		width *= rect[2]
		height *= rect[3]
		subax = fig.add_axes([x,y,width,height])
		
		x_labelsize = subax.get_xticklabels()[0].get_size()
		y_labelsize = subax.get_yticklabels()[0].get_size()
		x_labelsize *= rect[2]**0.5
		y_labelsize *= rect[3]**0.5
		
		subax.xaxis.set_tick_params(labelsize=x_labelsize)
		subax.yaxis.set_tick_params(labelsize=y_labelsize)
		
		return subax


	experimental_data = read_experimental_data_object(TAPobject_data.data_name)
	parameter_steps = pd.read_csv(TAPobject_data.output_name+'/optimization_results.csv')

	def generate_plot(k):
		#for k in range(0,parameter_steps.shape[0]):
		
		for j in parameter_steps.columns:
			if '[' in j:
				parser = copy.deepcopy(j)
				parse  = parser.split('[')[1]
				parse_num = parse.split(']')[0]

				if 'k' in parser:
					if 'f' in parser:
						TAPobject_data.mechanism.processes[int(parse_num)].f.k = parameter_steps[j][k]
						
					if 'b' in parser:
						TAPobject_data.mechanism.processes[int(parse_num)].b.k = parameter_steps[j][k]

		forward_problem(0.6,1,TAPobject_data)



		fig2, ax2 = plt.subplots(figsize=(10,5))
		##fig2 = Figure(figsize=(5, 4), dpi=100)
		##canvas = FigureCanvasAgg(fig2)
		##ax2 = fig2.gca()
		#fig2 = Figure()
		#canvas = fig2.canvas
		#ax2 = fig2.gca()

		ax2.set_xlabel('$Time\ (s)$', fontsize = 14)
	
		colors = ['b','orange','g','r','k','y','c','m','brown','darkgreen','goldenrod','lavender','lime']
	
		TAPobject_data.species.reference_diffusions = TAPobject_data.reactor.reference_diffusions
		zone_keys = list(TAPobject_data.species.reference_diffusions.keys())
	
		outlet_diffusions = {}
		for k in TAPobject_data.species.gasses.keys():
			outlet_diffusions[k] = TAPobject_data.species.reference_diffusions[zone_keys[-1]]*(mp.sqrt(TAPobject_data.species.reference_mass*TAPobject_data.species.temperature)/mp.sqrt(TAPobject_data.species.reference_temperature*TAPobject_data.species.gasses[k].mass))
		for k in TAPobject_data.species.inert_gasses.keys():
			outlet_diffusions[k] = TAPobject_data.species.reference_diffusions[zone_keys[-1]]*(mp.sqrt(TAPobject_data.species.reference_mass*TAPobject_data.species.temperature)/mp.sqrt(TAPobject_data.species.reference_temperature*TAPobject_data.species.inert_gasses[k].mass))
	
		if TAPobject_data.scaled_graph == True:
			ax2.set_ylabel('$Outlet\ Flow\ (1/s)$', fontsize = 14)
		else:
			ax2.set_ylabel('$Outlet\ Flow\ (nmol/s)$', fontsize = 14)
	
		try:
			#TAPobject_data.experimental_data = #output_data
			experimental_data = read_experimental_data_object(TAPobject_data.data_name)
			experimental_data_exist = True
		except:
			print('no experimental data')
			experimental_data_exist = False
		synthetic_data = read_experimental_data_object('./'+TAPobject_data.output_name+'/TAP_experimental_data.json')
		
		legend_label = []
		for j in TAPobject_data.species.gasses:
			legend_label.append(j)
		for j in TAPobject_data.species.inert_gasses:
			legend_label.append(j)
	
		for jnum,j in enumerate(TAPobject_data.species.gasses):
			if TAPobject_data.pulses_graphed == []:
				for k in synthetic_data[j]:
					if k == 0:
						plt.plot(synthetic_data['time'][0], synthetic_data[j][int(k)],label=j,color=colors[jnum],linestyle='--')
					else:
						plt.plot(synthetic_data['time'][0], synthetic_data[j][int(k)],color=colors[jnum],linestyle='--')
			else:
				for k in TAPobject_data.pulses_graphed:
					if k == 0:
						plt.plot(synthetic_data['time'][0], synthetic_data[j][int(k)],label=j,color=colors[jnum],linestyle='--')
					else:
						plt.plot(synthetic_data['time'][0], synthetic_data[j][int(k)],color=colors[jnum],linestyle='--')
				
	
		for jnum,j in enumerate(TAPobject_data.species.inert_gasses):
			if TAPobject_data.pulses_graphed == []:
				for k in synthetic_data[j]:
					if k == 0:
						plt.plot(synthetic_data['time'][0], synthetic_data[j][int(k)],label=j,color=colors[jnum+len(TAPobject_data.species.gasses.keys())],linestyle='--')
					else:
						plt.plot(synthetic_data['time'][0], synthetic_data[j][int(k)],color=colors[jnum+len(TAPobject_data.species.gasses.keys())],linestyle='--')
			else:
				for k in TAPobject_data.pulses_graphed:
					if k == 0:
						plt.plot(synthetic_data['time'][0], synthetic_data[j][int(k)],label=j,color=colors[jnum+len(TAPobject_data.species.gasses.keys())],linestyle='--')
					else:
						plt.plot(synthetic_data['time'][0], synthetic_data[j][int(k)],color=colors[jnum+len(TAPobject_data.species.gasses.keys())],linestyle='--')
				
		plt.ylim(0,4.6)
		plt.xlim(0,synthetic_data['time'][0][-1])
		if experimental_data_exist == True:
			for jnum,j in enumerate(TAPobject_data.species.gasses):#
				for k in experimental_data[j]:
					if k == 0:
						plt.scatter(experimental_data['time'][0][::6], experimental_data[j][0][::6],s=5,label=j+'_exp',color=colors[jnum])
					else:
						plt.scatter(experimental_data['time'][0][::6], experimental_data[j][0][::6],s=5,color=colors[jnum])
		
			for jnum,j in enumerate(TAPobject_data.species.inert_gasses):
				for k in experimental_data[j]:
					if k == 0:
						plt.scatter(experimental_data['time'][0][::6], experimental_data[j][0][::6],s=5,label=j+'_exp',color=colors[jnum+len(TAPobject_data.species.gasses.keys())])#
					else:
						plt.scatter(experimental_data['time'][0][::6], experimental_data[j][0][::6],s=5,color=colors[jnum+len(TAPobject_data.species.gasses.keys())])
	
		if TAPobject_data.display_analytical == True:
			# For reactant / product species
	
			analyticalTiming = synthetic_data['time'][0]#np.arange(0, time_steps*dt, dt).tolist()
			for jnum,kjc in enumerate(TAPobject_data.species.gasses):
				outlet = []
			
				factor = TAPobject_data.species.gasses[kjc].intensity*TAPobject_data.species.reference_pulse_size
				
				for k in analyticalTiming:
					analyticalValue = 0
					if k - TAPobject_data.species.gasses[kjc].delay > 0:
						for n in range(0,50):### Analytical solution section # gasses[name].catalyst_diffusion
							analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k - TAPobject_data.species.gasses[kjc].delay)*(outlet_diffusions[kjc]/(TAPobject_data.reactor.voids[0]*(TAPobject_data.reactor.length**2)))))
					else: 
						analyticalValue = -1
					if analyticalValue < 0:
						outlet.append(0)
					else:
						outlet.append(factor*(outlet_diffusions[kjc]*3.14159/(TAPobject_data.reactor.voids[0]*(TAPobject_data.reactor.length**2)))*analyticalValue)
				
				try:
					ax2.plot(analyticalTiming,outlet,color=colors[jnum],label='Analytical '+kjc, alpha=0.7)
				except ValueError:
					outlet = outlet[:-1]
					ax2.plot(analyticalTiming,outlet,color=colors[jnum],label='Analytical '+kjc, alpha=0.7)
	
			analyticalTiming = synthetic_data['time'][0]#np.arange(0, time_steps*dt, dt).tolist()
			for kjc in TAPobject_data.species.inert_gasses:
				outlet = []
			
				factor = TAPobject_data.species.inert_gasses[kjc].intensity*TAPobject_data.species.reference_pulse_size
				
				for k in analyticalTiming:
					analyticalValue = 0
					if k - TAPobject_data.species.inert_gasses[kjc].delay > 0:
						for n in range(0,50):### Analytical solution section # gasses[name].catalyst_diffusion
							analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k - TAPobject_data.species.inert_gasses[kjc].delay)*(outlet_diffusions[kjc]/(TAPobject_data.reactor.voids[0]*(TAPobject_data.reactor.length**2)))))
					else: 
						analyticalValue = -1
					if analyticalValue < 0:
						outlet.append(0)
					else:
						outlet.append(factor*(outlet_diffusions[kjc]*3.14159/(TAPobject_data.reactor.voids[0]*(TAPobject_data.reactor.length**2)))*analyticalValue)
				
				try:
					ax2.plot(analyticalTiming,outlet,color='k',label='Analytical '+kjc, alpha=0.7)
				except ValueError:
					outlet = outlet[:-1]
					ax2.plot(analyticalTiming,outlet,color='k',label='Analytical '+kjc, alpha=0.7)

		plt.legend()

		fig2.canvas.draw()

		#canvas.draw()
		buf = fig2.canvas.buffer_rgba()
		# convert to a NumPy array
		X = np.asarray(buf)
		#image = fig2
		#print(fig2.canvas.get_width_height()[::-1])
		#image  = image.reshape((int(6000000/4),2,2))#fig2.canvas.get_width_height()[::-1] + (6,))
	
		return X
		#	if TAPobject_data.mechanism.processes[k].f.use == 'G':
		#		TAPobject_data.mechanism.processes[k].f.Ga = Constant(TAPobject_data.mechanism.processes[k].f.Ga)
		#		TAPobject_data.mechanism.processes[k].f.dG = Constant(TAPobject_data.mechanism.processes[k].f.dG)
		#	elif TAPobject_data.mechanism.processes[k].f.use == 'E':
		#		TAPobject_data.mechanism.processes[k].f.Ao = Constant(TAPobject_data.mechanism.processes[k].f.Ao)
		#		TAPobject_data.mechanism.processes[k].f.Ea = Constant(TAPobject_data.mechanism.processes[k].f.Ea)
		#	elif TAPobject_data.mechanism.processes[k].f.use == 'k':
		#		TAPobject_data.mechanism.processes[k].f.k = Constant(TAPobject_data.mechanism.processes[k].f.k)		
		#	if TAPobject_data.mechanism.processes[k].b.use == 'E':
		#		TAPobject_data.mechanism.processes[k].b.Ao = Constant(TAPobject_data.mechanism.processes[k].b.Ao)
		#		TAPobject_data.mechanism.processes[k].b.Ea = Constant(TAPobject_data.mechanism.processes[k].b.Ea)
		#	elif TAPobject_data.mechanism.processes[k].b.use == 'k':
		#		try:
		#			TAPobject_data.mechanism.processes[k].b.k = Constant(TAPobject_data.mechanism.processes[k].b.k)
		#		except:
		#			pass

	kwargs_write = {'fps':4.0, 'quantizer':'nq'}
	components = list(range(0,parameter_steps.shape[0]))

	imageio.mimsave(TAPobject_data.output_name+'/output.gif', [generate_plot(i) for i in components], fps=4)
	#

	sys.exit()

def fluxGif(TAPobject_data: TAPobject):
	"""
	Return a gif showing changes made during the optimization process
	"""
	
	molecules = []
	
	for j in TAPobject_data.species.gasses.keys():
		molecules.append(j)
	try:
		experimental_data = read_experimental_data_object(TAPobject_data.data_name)
	except:
		pass

	colors = ['b','orange','g','r','k','y','c','m','brown','darkgreen','goldenrod','lavender','lime']
	
	TAPobject_data.species.reference_diffusions = TAPobject_data.reactor.reference_diffusions
	zone_keys = list(TAPobject_data.species.reference_diffusions.keys())
	
	outlet_diffusions = {}
	for k in TAPobject_data.species.gasses.keys():
		outlet_diffusions[k] = TAPobject_data.species.reference_diffusions[zone_keys[-1]]*(mp.sqrt(TAPobject_data.species.reference_mass*TAPobject_data.species.temperature)/mp.sqrt(TAPobject_data.species.reference_temperature*TAPobject_data.species.gasses[k].mass))
	for k in TAPobject_data.species.inert_gasses.keys():
		outlet_diffusions[k] = TAPobject_data.species.reference_diffusions[zone_keys[-1]]*(mp.sqrt(TAPobject_data.species.reference_mass*TAPobject_data.species.temperature)/mp.sqrt(TAPobject_data.species.reference_temperature*TAPobject_data.species.inert_gasses[k].mass))

	try:
		#TAPobject_data.experimental_data = #output_data
		experimental_data = read_experimental_data_object(TAPobject_data.data_name)
		experimental_data_exist = True
	except:
		print('no experimental data')
		experimental_data_exist = False
	synthetic_data = read_experimental_data_object('./'+TAPobject_data.output_name+'/TAP_experimental_data.json')
		
	legend_label = []
	for j in TAPobject_data.species.gasses:
		legend_label.append(j)
	for j in TAPobject_data.species.inert_gasses:
		legend_label.append(j)
		
	def generate_plot(kl):

		fig2, ax2 = plt.subplots(figsize=(10,5))
		ax2.set_xlabel('$Time\ (s)$', fontsize = 14)

		if TAPobject_data.scaled_graph == True:
			ax2.set_ylabel('$Outlet\ Flow\ (1/s)$', fontsize = 14)
		else:
			ax2.set_ylabel('$Outlet\ Flow\ (nmol/s)$', fontsize = 14)

		for jnum,j in enumerate(TAPobject_data.species.gasses):
			if TAPobject_data.pulses_graphed == []:
				plt.plot(synthetic_data['time'][0], synthetic_data[j][kl],label=j,color=colors[jnum],linestyle='--')	
	
		for jnum,j in enumerate(TAPobject_data.species.inert_gasses):
			if TAPobject_data.pulses_graphed == []:
				plt.plot(synthetic_data['time'][0], synthetic_data[j][kl],label=j,color=colors[jnum+len(TAPobject_data.species.gasses.keys())],linestyle='--')
				
		plt.xlim(0,1)
		plt.ylim(0,6)
		if experimental_data_exist == True:
			for jnum,j in enumerate(TAPobject_data.species.gasses):#
				plt.scatter(experimental_data['time'][0][::6], experimental_data[j][kl][::6],s=5,label=j+'_exp',color=colors[jnum])
		
			for jnum,j in enumerate(TAPobject_data.species.inert_gasses):
				plt.scatter(experimental_data['time'][0][::6], experimental_data[j][kl][::6],s=5,label=j+'_exp',color=colors[jnum+len(TAPobject_data.species.gasses.keys())])#
				
		if TAPobject_data.display_analytical == True:
			# For reactant / product species
	
			analyticalTiming = synthetic_data['time'][0]#np.arange(0, time_steps*dt, dt).tolist()
			for jnum,kjc in enumerate(TAPobject_data.species.gasses):
				outlet = []
			
				factor = TAPobject_data.species.gasses[kjc].intensity*TAPobject_data.species.reference_pulse_size
				
				for k in analyticalTiming:
					analyticalValue = 0
					if k - TAPobject_data.species.gasses[kjc].delay > 0:
						for n in range(0,50):### Analytical solution section # gasses[name].catalyst_diffusion
							analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k - TAPobject_data.species.gasses[kjc].delay)*(outlet_diffusions[kjc]/(TAPobject_data.reactor.voids[0]*(TAPobject_data.reactor.length**2)))))
					else: 
						analyticalValue = -1
					if analyticalValue < 0:
						outlet.append(0)
					else:
						outlet.append(factor*(outlet_diffusions[kjc]*3.14159/(TAPobject_data.reactor.voids[0]*(TAPobject_data.reactor.length**2)))*analyticalValue)
				
				try:
					ax2.plot(analyticalTiming,outlet,color=colors[jnum],label='Analytical '+kjc, alpha=0.7)
				except ValueError:
					outlet = outlet[:-1]
					ax2.plot(analyticalTiming,outlet,color=colors[jnum],label='Analytical '+kjc, alpha=0.7)
	
			analyticalTiming = synthetic_data['time'][0]#np.arange(0, time_steps*dt, dt).tolist()
			for kjc in TAPobject_data.species.inert_gasses:
				outlet = []
			
				factor = TAPobject_data.species.inert_gasses[kjc].intensity*TAPobject_data.species.reference_pulse_size
				
				for k in analyticalTiming:
					analyticalValue = 0
					if k - TAPobject_data.species.inert_gasses[kjc].delay > 0:
						for n in range(0,50):### Analytical solution section # gasses[name].catalyst_diffusion
							analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k - TAPobject_data.species.inert_gasses[kjc].delay)*(outlet_diffusions[kjc]/(TAPobject_data.reactor.voids[0]*(TAPobject_data.reactor.length**2)))))
					else: 
						analyticalValue = -1
					if analyticalValue < 0:
						outlet.append(0)
					else:
						outlet.append(factor*(outlet_diffusions[kjc]*3.14159/(TAPobject_data.reactor.voids[0]*(TAPobject_data.reactor.length**2)))*analyticalValue)
				
				try:
					ax2.plot(analyticalTiming,outlet,color='k',label='Analytical '+kjc, alpha=0.7)
				except ValueError:
					outlet = outlet[:-1]
					ax2.plot(analyticalTiming,outlet,color='k',label='Analytical '+kjc, alpha=0.7)

		plt.legend()

		fig2.canvas.draw()

		#canvas.draw()
		buf = fig2.canvas.buffer_rgba()
		# convert to a NumPy array
		X = np.asarray(buf)
		#image = fig2
		#print(fig2.canvas.get_width_height()[::-1])
		#image  = image.reshape((int(6000000/4),2,2))#fig2.canvas.get_width_height()[::-1] + (6,))
	
		return X

	pulses_observed = []
	for k in synthetic_data[j]:
		pulses_observed.append(int(k))


	kwargs_write = {'fps':4.0, 'quantizer':'nq'}
	pulses_observed.sort()#list(range(0,len(pulses_observed)))


	imageio.mimsave(TAPobject_data.output_name+'/flux_visual.gif', [generate_plot(i) for i in pulses_observed], fps=8)
	#

	sys.exit()

def concentrationGif(TAPobject_data: TAPobject,pulse_number=None):
	"""
	Return a gif showing changes made during the optimization process
	"""
	
	colors = ['b','orange','g','r','k','y','c','m','brown','darkgreen','goldenrod','lavender','lime']
	
	for knum,k in enumerate(TAPobject_data.conc_profiles):
		def generate_plot(lnum):
			
			mesh_data = pd.read_csv(TAPobject_data.output_name+'/mesh.csv',header=None)
			experimental_data = pd.read_csv(TAPobject_data.output_name+'/'+k+'_pulse_'+str(lnum)+'.csv',header=None)
			fig2, ax2 = plt.subplots(figsize=(10,5))
			ax2.set_xlabel('$Reactor Length\ (cm)$', fontsize = 14)
			ax2.set_ylabel('$Concentration\ (nmol/cm3)$', fontsize = 14)
			ax2.set_ylim(0,10)
			plt.plot(mesh_data[0][:].transpose(), experimental_data.iloc[1][:],label=knum,color=colors[knum])	
			
			fig2.canvas.draw()

			#canvas.draw()
			buf = fig2.canvas.buffer_rgba()
			# convert to a NumPy array
			X = np.asarray(buf)
			#image = fig2
			#print(fig2.canvas.get_width_height()[::-1])
			#image  = image.reshape((int(6000000/4),2,2))#fig2.canvas.get_width_height()[::-1] + (6,))
			
			return X
			
		pulses_observed = []
		for lnum in range(0,pulse_number):
			pulses_observed.append(int(lnum))

		kwargs_write = {'fps':.0, 'quantizer':'nq'}
		pulses_observed.sort()#list(range(0,len(pulses_observed)))

		imageio.mimsave(TAPobject_data.output_name+'/surface_'+k+'.gif', [generate_plot(i) for i in pulses_observed[::10]], fps=8)
	#
		sys.exit()

def single_pulse_concentrationGif(TAPobject_data: TAPobject,pulse_number=None):
	"""
	Return a gif showing changes made during the optimization process
	"""
	
	colors = ['b','orange','g','r','k','y','c','m','brown','darkgreen','goldenrod','lavender','lime']
	
	for knum,k in enumerate(TAPobject_data.conc_profiles):
		mesh_data = pd.read_csv(TAPobject_data.output_name+'/mesh.csv',header=None)
		experimental_data = pd.read_csv(TAPobject_data.output_name+'/'+k+'_pulse_'+str(pulse_number)+'.csv',header=None)
			
		def generate_plot(lnum):
			
			fig2, ax2 = plt.subplots(figsize=(10,5))
			ax2.set_xlabel('$Reactor Length\ (cm)$', fontsize = 14)
			ax2.set_ylabel('$Concentration\ (nmol/cm3)$', fontsize = 14)
			if k == '*(0)':
				ax2.set_ylim(0,10)
			elif k == 'C':
				ax2.set_ylim(0,0.5)
			elif k == 'D*(0)':
				ax2.set_ylim(0,2)
			else:
				ax2.set_ylim(0,50)
			
			plt.plot(mesh_data[0][:].transpose(), experimental_data.iloc[lnum][:],label=knum,color=colors[knum])	
			
			fig2.canvas.draw()

			#canvas.draw()
			buf = fig2.canvas.buffer_rgba()
			# convert to a NumPy array
			X = np.asarray(buf)
			#image = fig2
			#print(fig2.canvas.get_width_height()[::-1])
			#image  = image.reshape((int(6000000/4),2,2))#fig2.canvas.get_width_height()[::-1] + (6,))
			
			return X

		pulses_observed = []
		for lnum in range(0,1500):
			pulses_observed.append(int(lnum))

		kwargs_write = {'fps':34.0, 'quantizer':'nq'}
		pulses_observed.sort()#list(range(0,len(pulses_observed)))

		imageio.mimsave(TAPobject_data.output_name+'/individual_profile_'+k+'_'+str(pulse_number)+'.gif', [generate_plot(i) for i in pulses_observed], fps=8)
	#
	sys.exit()

	
	