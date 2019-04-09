import pandas as pd
import matplotlib.pyplot as plt
import sys

data = pd.ExcelFile('./COoxidationData_Pt/COoxidationData.xlsx')

df1 = pd.read_excel(data,'CO',header=None)
df2 = pd.read_excel(data,'CO2',header=None)
df3 = pd.read_excel(data,'Argon',header=None)

#ax2.plot(df1[0].convert_objects(convert_numeric=True),df1[1].convert_objects(convert_numeric=True))
#ax2.plot(df2[0].convert_objects(convert_numeric=True),df2[1].convert_objects(convert_numeric=True))
#ax2.plot(df3[0].convert_objects(convert_numeric=True),df3[1].convert_objects(convert_numeric=True))

#plt.ion()

#plt.show()
#x = df1[0].convert_objects(convert_numeric=True)
#y = df1[0].convert_objects(convert_numeric=True)
#line1, = ax2.plot(x, y, 'b-')
#line1, = ax.plot(x, y, 'b-')

for k in range(1,50):
	fig2,ax2 = plt.subplots()
	#ax2.plot(df1[0].convert_objects(convert_numeric=True),df1[k].convert_objects(convert_numeric=True))
	#ax2.plot(df2[0].convert_objects(convert_numeric=True),df2[k+1].convert_objects(convert_numeric=True))
	ax2.plot(df3[0][:2000].convert_objects(convert_numeric=True),df3[k][:2000].convert_objects(convert_numeric=True))
	#line1.set_ydata(df1[k].convert_objects(convert_numeric=True))
	plt.show()
#fig2.canvas.draw()


	#plt.clf()
	#fig2.clear()
	
#plt.show()

#ax2.plot(df1[0].convert_objects(convert_numeric=True),df1[50].convert_objects(convert_numeric=True))
#ax2.plot(df2[0].convert_objects(convert_numeric=True),df2[50].convert_objects(convert_numeric=True))
#ax2.plot(df3[0].convert_objects(convert_numeric=True),df3[50].convert_objects(convert_numeric=True))

#plt.show()