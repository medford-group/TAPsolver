def molecularProperties(gasSpecies,propValue,temperature=398):

	def thermoCalculation(T,dH,A,B,C,D,E,F,G,H):
		temp = T/1000

		h = dH + (A*temp + (B*(temp**2)/2) + (C*(temp**3)/3) + (D*(temp**4)/4) - (E/temp) + F - H)#/1000
		s = A*math.log(temp) + B*temp + C*(temp**2)/2 + D*(temp**3)/3 - E/(2*temp**2) + G

		return h,s/1000

	molecularValues = {}

	## Nist Chemistry WebBook was used to determine the shomate equation parameters
	## 

	# Oxygen Scrambling
	molecularValues['O218'] = {'mass':36,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	molecularValues['O216'] = {'mass':32,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	molecularValues['O2'] = {'mass':32,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	molecularValues['O18O16'] = {'mass':34,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	molecularValues['O16O18'] = {'mass':34,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	
	# Natural Gas
	molecularValues['CO'] = {'mass':28.01,'shomate':{'dH':-110.5271,'A':25.5679,'B':6.096130,'C':4.05465,'D':-2.67130,'E':0.131021,'F':-118.0089,'G':227.366,'H':-110.5271,'Trange':[298,1300]}}
	molecularValues['CO2'] = {'mass':44.01,'shomate':{'dH':-393.5224,'A':24.99735,'B':55.18696,'C':-33.69137,'D':7.948387,'E':-0.136638,'F':-403.6075,'G':228.2431,'H':-393.5224,'Trange':[298,1200]}}
	#molecularValues['CH2'] = {'mass':14.01}
	#molecularValues['CH3'] = {'mass':15.01}
	molecularValues['CH4'] = {'mass':16.01,'shomate':{'dH':-74.87310,'A':-0.703029,'B':108.4773,'C':-42.52157,'D':5.862788,'E':0.678565,'F':-76.84376,'G':158.7163,'H':-74.87310,'Trange':[298,1300]}}
	molecularValues['C2H4'] = {'mass':32,'shomate':{'dH':52.46694,'A':-6.387880,'B':184.4019,'C':-112.9718,'D':28.49593,'E':0.315540,'F':48.17332,'G':163.1568,'H':52.46694,'Trange':[298,1200]}}
	#molecularValues['C2H6'] = {'mass':34,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C3H4'] = {'mass':42.03,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C3H6'] = {'mass':42.03,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C3H8'] = {'mass':44.03,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C4H6'] = {'mass':54.04,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C4H8'] = {'mass':56.04,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C4H10'] = {'mass':58.04,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C5H10'] = {'mass':70.05,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C5H12'] = {'mass':72.05,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	molecularValues['H2S'] = {'mass':34,'shomate':{'dH':-20.502,'A':26.88412	,'B':18.67809,'C':3.434203,'D':-3.378702,'E':0.135882	,'F':-28.91211,'G':233.3747,'H':-20.502,'Trange':[298,1400]}}
	#molecularValues['H2O'] = {'mass':18,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	molecularValues['SO2'] = {'mass':64,'shomate':{'dH':-296.8422,'A':21.43049,'B':74.35094,'C':-57.75217,'D':16.35534,'E':0.086731,'F':-305.7688,'G':254.8872,'H':-296.8422,'Trange':[298,1200]}}

	# Atmospheric / Exhaust / Fertilizer
	molecularValues['N2'] = {'mass':28,'shomate':{'dH':0.0,'A':28.98641,'B':1.853978,'C':-9.647459,'D':16.63537,'E':0.000117,'F':-8.671914,'G':226.4168,'H':0.0,'Trange':[100,500]}}
	molecularValues['NO'] = {'mass':30,'shomate':{'dH':90.29114,'A':23.83491,'B':12.58878,'C':-1.139011,'D':-1.497459,'E':0.214194,'F':83.35783,'G':237.1219,'H':90.29114,'Trange':[298,1200]}}
	molecularValues['NO2'] = {'mass':46,'shomate':{'dH':33.09502,'A':16.10857,'B':75.89525,'C':-54.38740,'D':14.30777,'E':0.239423,'F':26.17464,'G':240.5386,'H':33.09502,'Trange':[298,1200]}}
	molecularValues['N2O'] = {'mass':44.01,'shomate':{'dH':82.04824,'A':27.67988,'B':51.14898,'C':-30.64454,'D':6.847911,'E':-0.157906,'F':71.24934,'G':238.6164,'H':82.04824,'Trange':[298,1400]}}
	molecularValues['NO3'] = {'mass':62,'shomate':{'dH':71.12800,'A':11.22316	,'B':166.3889,'C':-148.4458,'D':47.40598,'E':-0.176791,'F':61.00858	,'G':221.7679,'H':71.12800,'Trange':[298,1100]}}
	molecularValues['N2O5'] = {'mass':108,'shomate':{'dH':82.84320,'A':39.09663,'B':114.8006,'C':-81.97125,'D':21.77249,'E':-0.088738,'F':66.46786,'G':324.5776,'H':82.84320,'Trange':[298,1100]}}
	molecularValues['NH3'] = {'mass':17.03,'shomate':{'dH':-45.89806,'A':19.99563,'B':49.77119,'C':-15.37599,'D':1.921168,'E':0.189174,'F':-53.30667,'G':203.8591,'H':-45.89806,'Trange':[298,1400]}}
	molecularValues['HNO3'] = {'mass':63,'shomate':{'dH':-134.3060,'A':19.63229,'B':153.9599,'C':-115.8378,'D':32.87955,'E':-0.249114,'F':-146.8818,'G':247.7049,'H':-134.3060,'Trange':[298,1200]}}
	molecularValues['O3'] = {'mass':48,'shomate':{'dH':142.6740,'A':21.66157,'B':79.86001,'C':-66.02603,'D':19.58363,'E':-0.079251,'F':132.9407,'G':243.6406,'H':142.6740,'Trange':[298,1200]}}
	molecularValues['N2H4'] = {'mass':32}


	# Other diatomic / common species
	molecularValues['Cl2'] = {'mass':71,'shomate':{'dH':0.0,'A':33.05060,'B':12.22940,'C':-12.06510,'D':4.385330,'E':-0.159494,'F':-10.83480,'G':259.0290,'H':0.0,'Trange':[298,1000]}}
	molecularValues['Br2'] = {'mass':159.8,'shomate':{'dH':30.91001,'A':38.52723,'B':-1.976835,'C':1.526107,'D':-0.198398,'E':-0.185815,'F':18.87620,'G':291.4863,'H':30.91001,'Trange':[333,3400]}}
	molecularValues['I2'] = {'mass':253.8,'shomate':{'dH':62.42110,'A':37.79763	,'B':0.225453,'C':-0.912556,'D':1.034913,'E':-0.083826,'F':50.86865,'G':305.9199,'H':62.42110,'Trange':[458,2000]}}
	molecularValues['H2'] = {'mass':2.01,'shomate':{'dH':0.0,'A':33.066178,'B':-11.363417,'C':11.432816,'D':-2.772874,'E':-0.158558,'F':-9.980797,'G':172.707974,'H':0.0,'Trange':[298,1000]}}
	molecularValues['HCl'] = {'mass':36.5,'shomate':{'dH':-92.312,'A':32.12392,'B':-13.45805,'C':19.86852,'D':-6.853936,'E':-0.049672,'F':-101.6206,'G':228.6866,'H':-92.312,'Trange':[298,1200]}}

	molecularValues['He'] = {'mass':4,'shomate':{'dH':0.0,'A':20.78603,'B':4.8506e-10,'C':-1.5829e-10,'D':1.5251e-11,'E':3.1963e-11,'F':-6.197341,'G':151.3064,'H':0.0,'Trange':[298,6000]}}
	molecularValues['Ne'] = {'mass':20.18,'shomate':{'dH':0.0,'A':20.78603,'B':4.8506e-10,'C':-1.5829e-10,'D':1.525102e-11,'E':3.1963e-11,'F':-6.197341,'G':171.48,'H':0.0,'Trange':[298,6000]}}
	molecularValues['Ar'] = {'mass':40,'shomate':{'dH':0.0,'A':20.786,'B':2.825911e-7,'C':-1.46419e-7,'D':1.092131e-8,'E':-3.661371e-8,'F':-6.197350,'G':179.999,'H':0.0,'Trange':[298,6000]}}
	molecularValues['Kr'] = {'mass':83.798,'shomate':{'dH':0.0,'A':20.78603,'B':4.850638e-10,'C':-1.582916e-10,'D':1.525102e-11,'E':3.196347e-11,'F':-6.197341,'G':189.239,'H':0.0,'Trange':[298,6000]}}
	molecularValues['Xe'] = {'mass':131.293,'shomate':{'dH':0.0,'A':20.786,'B':7.4493e-7,'C':-2.0494e-7,'D':1.066e-8,'E':2.500261e-8,'F':-6.197350,'G':194.8380,'H':0.0,'Trange':[298,6000]}}

	if gasSpecies in molecularValues.keys():
		if propValue == 'mass':
			return molecularValues[gasSpecies][propValue]
		elif propValue == 'freeEnergy':

			tempShomate = molecularValues[gasSpecies]['shomate']
			hnew,snew = thermoCalculation(temperature,tempShomate['dH'],tempShomate['A'],tempShomate['B'],tempShomate['C'],tempShomate['D'],tempShomate['E'],tempShomate['F'],tempShomate['G'],tempShomate['H'])
			
			return hnew-temperature*snew
	else:
		return ''