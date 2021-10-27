def progressBar(value, endvalue, bar_length=20):
	
	""" Generate the progress bar"""

	percent = float(value) / endvalue
	arrow = '-' * int(round(percent * bar_length)-1) + '>'
	spaces = ' ' * (bar_length - len(arrow))
	sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, round(percent * 100,3)))
	sys.stdout.flush()

def processTime(start_time):
	if (time.time() - start_time) < 120:
		return 'Completed in: '+str(round((time.time() - start_time),3))+' seconds'
	elif (time.time() - start_time)/60 < 120:
		return 'Completed in: '+str(round((time.time() - start_time)/60,3))+' minutes'
	else:
		return 'Completed in: '+str(round((time.time() - start_time)/3600,3))+' hours'
