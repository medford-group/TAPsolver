def establishOutletGraph(reactor,gas_phase,reacs,inerts,scaleGraph):

	"""Generate the outlet flux graph (initialize, anyway)"""

	fig2, ax2 = plt.subplots()
	ax2.set_xlabel('$Time\ (s)$', fontsize = 14)

	if reactor == 'tap':
		if scaleGraph.lower() == 'true':
			ax2.set_ylabel('$Outlet\ Flow\ (1/s)$', fontsize = 14)
		else:
			ax2.set_ylabel('$Outlet\ Flow\ (nmol/s)$', fontsize = 14)
	elif reactor == 't_pfr' or 't_pfr_diff':
		ax2.set_ylabel('$Outlet Concentration (molecules/cm3)$')

	elif reactor == 'batch':
		print('batch')
		ax2.set_ylabel('$Concentration$')

	legend_label = []
	header = "time"
	for k in range(0,gas_phase):
		legend_label.append(reacs[k])
		header = header+","+reacs[k]
	if reactor != 'batch':
		for j in range(0,inerts):
			legend_label.append("Inert-"+str(1+j))
		header = header+",Inert"
	
	colors = ['b','orange','g','r','k','y','c','m','brown','darkgreen','goldenrod','lavender','lime']

	return fig2,ax2,legend_label,header,colors