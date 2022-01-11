import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats
import math
import imageio
import matplotlib.pylab as pl

time_steps = 2000
dt = 2/2000
colors = pl.cm.jet(np.linspace(0,1,30))
analyticalTiming = np.arange(0, time_steps*dt, dt).tolist()
outlet = []
factor = float(1)*1

outlet_3 = [0]

for j in analyticalTiming:
	#outlet_3.append(10*np.log(1+0.02*j))
	outlet_3.append(0.3+(j/2)/(1.5-(j/2)))
outlet_3 = np.array(outlet_3)	
outlet_3 = np.flip(outlet_3,axis=0)		
for k in range(0,2000+1):
	analyticalValue = 0.0
	if k*dt > 0:
		for n in range(0,50):
			analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt)*(16/(0.4*(5**2)))))
	else: 
		analyticalValue = -1
	if analyticalValue < 0:
		outlet.append(0)
	else:
		outlet.append(factor*(16*3.14159/(0.4*(5**2)))*analyticalValue)
#	try:
#		ax2.plot(analyticalTiming,outlet,color='k', alpha=(1-step)/total_steps)
#	except ValueError:
#		outlet = outlet[:-1]
#		ax2.plot(analyticalTiming,outlet,color='k', alpha=(1-step)/total_steps)
outlet_2 = 0.6*stats.norm.pdf(analyticalTiming, 1, 0.15)					
outlet = np.array(outlet)
total_steps = 30
new_font = 'Georgia'
new_font_2 = 'Lucida Console'
def make_gif():
		
	"""
	Generates the subplot to help visualize some other quantitiy
	"""

	def tap_plot(step):
		fig, ax = plt.subplots(figsize=(10,6))
		hfont = {'fontname':'Helvetica'}
		#plt.rcParams["font.family"] = "cursive"
		ax.grid()
		#plt.tight_layout()
		props = dict(facecolor='white')
		
		if step < 30:
			ax.plot(analyticalTiming,outlet[:-1],color=colors[0],linewidth=6)
			ub_ci = []
			lb_ci = []
			for j_top in outlet[:-1]:
				ub_ci.append(j_top+0.5)
				lb_ci.append(j_top-0.5)
			ub_ci = np.array(ub_ci)
			lb_ci = np.array(lb_ci)
			ax.fill_between(analyticalTiming, lb_ci, ub_ci, color=colors[0], alpha=.35)
			ax.text(1.02, 3.35,'1. Run a TAP Experiment', fontsize=20,verticalalignment='top',alpha=1,bbox=props,fontname=new_font)
			ax.text(0.88, -0.32,'Time (s)', fontsize=20,verticalalignment='top',alpha=1,fontname=new_font_2)
			ax.text(-0.22, 2.4,'Flow (nmol/s)', rotation=90,fontsize=20,verticalalignment='top',alpha=1,fontname=new_font_2)
		elif 30 <= step and step <= 60:
			ax.text(1.02, 3.35,'1. Run a TAP Experiment', fontsize=20,verticalalignment='top',alpha=(60-step)/30,fontname=new_font)
			ax.text(1.02, 2.85,'2. Fit Parameters + UQ', fontsize=20,verticalalignment='top',alpha=(step-30)/30,fontname=new_font)
			
			ax.text(0.88, -0.32,'Time (s)', fontsize=20,verticalalignment='top',alpha=(60-step)/30,fontname=new_font_2)
			ax.text(0.65, -0.32,'Kinetic Parameter (1/s)', fontsize=20,verticalalignment='top',alpha=(step-30)/30,fontname=new_font_2)

			ax.text(-0.22, 2.4,'Flow (nmol/s)', rotation=90,fontsize=20,verticalalignment='top',alpha=(60-step)/30,fontname=new_font_2)
			ax.text(-0.22, 2.7,'Probability Density', rotation=90,fontsize=20,verticalalignment='top',alpha=(step-30)/30,fontname=new_font_2)

			output_temp = np.roll(outlet[:-1]*((30-(step-31))/30),shift=(step-30)*24)+outlet_2*((step-31)/30)
			ub_ci = []
			lb_ci = []
			for j_top in output_temp:
				ub_ci.append(j_top+0.5*((30-(step-31))/30))
				lb_ci.append(j_top-0.5*((30-(step-31))/30))
			ub_ci = np.array(ub_ci)
			lb_ci = np.array(lb_ci)
			ax.fill_between(analyticalTiming, lb_ci, ub_ci,color=colors[step-31], alpha=.35*(30-(step-30))/30)
			ax.plot(analyticalTiming,output_temp,color=colors[step-31],linewidth=5)
		elif 60 <= step and step <= 90:
			ax.text(1.02, 2.85,'2. Fit Parameters + UQ', fontsize=20,verticalalignment='top',bbox=props,fontname=new_font)
			ax.text(0.65, -0.32,'Kinetic Parameter (1/s)', fontsize=20,verticalalignment='top',alpha=1,fontname=new_font_2)
			ax.text(-0.22, 2.7,'Probability Density', rotation=90,fontsize=20,verticalalignment='top',alpha=1,fontname=new_font_2)
			ax.plot(analyticalTiming,outlet_2,color=colors[29],linewidth=6)
		elif 90 <= step and step <= 120:
			ax.text(1.02, 2.85,'2. Fit Parameters + UQ', fontsize=20,verticalalignment='top',alpha=(120-step)/30,fontname=new_font)
			ax.text(1.02, 2.35,'3. Propagate to Reactor', fontsize=20,verticalalignment='top',alpha=(step-90)/30,fontname=new_font)

			ax.text(0.65, -0.32,'Kinetic Parameter (1/s)', fontsize=20,verticalalignment='top',alpha=(120-step)/30,fontname=new_font_2)
			ax.text(-0.22, 2.7,'Probability Density', rotation=90,fontsize=20,verticalalignment='top',alpha=(120-step)/30,fontname=new_font_2)

			ax.text(0.72, -0.32,'Reactor Length (m)', fontsize=20,verticalalignment='top',alpha=(step-90)/30,fontname=new_font_2)
			ax.text(-0.22, 3,'Concentration (mol/cm3)', rotation=90,fontsize=20,verticalalignment='top',alpha=(step-90)/30,fontname=new_font_2)

			output_temp = np.roll(outlet_2*((30-(step-91))/30),shift=(91-step)*24)+outlet_3[:-1]*((step-91)/30)
			ub_ci = []
			lb_ci = []
			for j_top in output_temp:
				ub_ci.append(j_top+0.5*((step-91)/30))
				lb_ci.append(j_top-0.5*((step-91)/30))
			ub_ci = np.array(ub_ci)
			lb_ci = np.array(lb_ci)
			ax.fill_between(analyticalTiming, lb_ci, ub_ci, color=colors[120-step], alpha=.35*((step-91)/30))
			ax.plot(analyticalTiming,output_temp,color=colors[120-step],linewidth=6)

		
		elif 120 <= step and step <= 150:
			ax.text(1.02, 2.35,'3. Propagate to Reactor', fontsize=20,verticalalignment='top',bbox=props,fontname=new_font)
			ax.text(0.72, -0.32,'Reactor Length (m)', fontsize=20,verticalalignment='top',alpha=1,fontname=new_font_2)
			ax.text(-0.22, 3,'Concentration (mol/cm3)', rotation=90,fontsize=20,verticalalignment='top',alpha=1,fontname=new_font_2)
			ub_ci = []
			lb_ci = []
			for j_top in outlet_3[:-1]:
				ub_ci.append(j_top+0.5)
				lb_ci.append(j_top-0.5)
			ub_ci = np.array(ub_ci)
			lb_ci = np.array(lb_ci)
			ax.fill_between(analyticalTiming, lb_ci, ub_ci, color=colors[0], alpha=.35)
			ax.plot(analyticalTiming,outlet_3[:-1],color=colors[0],linewidth=6)

		elif 150 <= step and step <= 180:
			ax.text(1.02, 2.35,'3. Propagate to Reactor', fontsize=20,verticalalignment='top',alpha=(180-step)/30,fontname=new_font)
			ax.text(1.02, 3.35,'1. Run a TAP Experiment', fontsize=20,verticalalignment='top',alpha=(step-150)/30,fontname=new_font)

			ax.text(0.72, -0.32,'Reactor Length (m)', fontsize=20,verticalalignment='top',alpha=(180-step)/30,fontname=new_font_2)
			ax.text(-0.22, 3,'Concentration (mol/cm3)', rotation=90,fontsize=20,verticalalignment='top',alpha=(180-step)/30,fontname=new_font_2)

			ax.text(0.88, -0.32,'Time (s)', fontsize=20,verticalalignment='top',alpha=(step-150)/30,fontname=new_font_2)
			ax.text(-0.22, 2.4,'Flow (nmol/s)', rotation=90,fontsize=20,verticalalignment='top',alpha=(step-150)/30,fontname=new_font_2)

			output_temp = outlet_3[:-1]*((30-(step-151))/30)+outlet[:-1]*((step-151)/30)
			ub_ci = []
			lb_ci = []
			for j_top in output_temp:
				ub_ci.append(j_top+0.5)
				lb_ci.append(j_top-0.5)
			ub_ci = np.array(ub_ci)
			lb_ci = np.array(lb_ci)
			ax.fill_between(analyticalTiming, lb_ci, ub_ci, color=colors[0], alpha=.35)

			#ci = 1.96 * np.std(y)/np.sqrt(len(x))
			ax.plot(analyticalTiming,output_temp,color=colors[0],linewidth=6)

		#ax.set_xlabel('Time (s)', fontsize=16)
		ax.set_ylim(0,3.5)
		ax.set_xlim(0,2)
		plt.xticks(fontsize=20)
		plt.yticks(fontsize=20)
		#ax.set_ylabel('Flow (nmol/s)', fontsize=16)
		fig.canvas.draw()
		image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
		image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
		return image
	
	kwargs_write = {'fps':4.0, 'quantizer':'nq'}
	components = list(range(0,180))

	imageio.mimsave('./output.gif', [tap_plot(i) for i in components], fps=20)

make_gif()